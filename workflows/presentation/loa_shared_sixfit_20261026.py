"""Observation-only six-fit shared-state engineering tests and resource benchmark.

All source predictions remain untouched. No training, truth, tides, or production
publication. A benchmark at NFE2 is never a posterior product. Production needs
the resulting resource measurements and a separately frozen launch contract.
"""
import argparse
from dataclasses import replace
import gc
import hashlib
import itertools
import json
import os
from pathlib import Path
import socket
import sys
import time

os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
ILLUSTRIS = Path(os.environ.get('ILLUSTRIS_ROOT', '/global/u2/d/dkololgi/TNG/Illustris'))
sys.path.insert(0, str(ILLUSTRIS))
import numpy as np
import torch
from workflows.sbi import e2e_coupled_views as views
from workflows.sbi import e2e_cfm_pilot_evaluate as evaluate
from workflows.sbi import e2e_loa_field_sampling as settings
from workflows.sbi.e2e_global_canvas import Geometry, integrate
from workflows.sbi.e2e_coupled_benchmark_models import FieldCondition
from workflows.sbi.e2e_vdm_context_models import project, replicate, decode_density, block_mean
from workflows.sbi.e2e_wholefield_pilot_models import FineModel, Condition, sample
from workflows.sbi.e2e_wholefield_pilot_train import sources as training_sources

INPUT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1')
FAMILIES = ('legacy', 'unet', 'transformer')

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(8<<20),b''):h.update(b)
    return h.hexdigest()

def read(path):return json.loads(Path(path).read_text())

def write(path, value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False)

def geometry(midpoints, expanded=False):
    if not expanded:return Geometry.make(midpoints,'fine')
    mid=np.asarray(midpoints,dtype=int)
    if mid.ndim!=2 or mid.shape[1]!=3 or not len(mid) or np.any(mid%8):raise ValueError('Invalid raw-cell lattice')
    if len(np.unique(mid,axis=0))!=len(mid):raise ValueError('Duplicate windows')
    window=np.array([96,64,64]);starts=(mid-window)//2
    origin=starts.min(0);relative=starts-origin;shape=(relative+window).max(0)
    if np.any(relative%4) or np.any(shape%4):raise ValueError('Invalid shared block lattice')
    return Geometry(tuple(map(tuple,relative)),tuple(shape),tuple(origin),tuple(window))

class ExpandedAdapter:
    """Keep the trained endpoint/velocity conversion and BF16 forward unchanged."""
    def __init__(self, model):self.model=model
    def schedule(self,t):return t
    def __call__(self,x,t,condition):
        with torch.autocast('cuda',dtype=torch.bfloat16):return self.model(x,t,condition)

@torch.no_grad()
def integrate_expanded(model,conditions,noise,g,nfe=128,on_step=None,cache_conditions=True):
    """Native expanded-model Heun arithmetic with synchronous overlap consensus.

    Window velocities already have zero block mean. Coverage multiplicity is
    constant within each aligned block, so their weighted mean needs no second
    projection. An extra projection changes BF16 trajectories and fails native
    single-window parity. Preserve the native initial/proposal/update projections.
    """
    if nfe<2 or nfe%2 or len(conditions)!=len(g.starts) or tuple(noise.shape[2:])!=g.shape:
        raise ValueError('Invalid shared expanded integration contract')
    if any(n%4 for n in (*g.shape,*g.window)) or any(v%4 for p in g.starts for v in p):
        raise ValueError('Consensus requires complete aligned blocks')
    weights=g.weights(noise.device);denominator=weights.clamp_min(1)
    state=project(noise.float());model.eval()
    with torch.autocast('cuda',dtype=torch.bfloat16,enabled=state.is_cuda):
        # Small engineering panels fit cached conditions; full caps prepare one
        # window at a time to keep residency bounded on the approved GPU.
        prepared=[model.prepare(c) for c in conditions] if cache_conditions and len(conditions)<=8 else None
        def velocity(value,t):
            total=torch.zeros_like(value)
            for i,condition in enumerate(conditions):
                cut=g.crop(i);x=value[cut]
                p=prepared[i] if prepared is not None else model.prepare(condition)
                total[cut]+=model.velocity(x,x.new_full((len(x),),t),p)
            return total/denominator
        for step in range(nfe//2):
            dt=2/nfe;t=step*dt
            a=velocity(state,t);b=velocity(project(state+dt*a),t+dt)
            state=project(state+.5*dt*(a+b))
            if on_step is not None:on_step(step+1,state)
    if not torch.isfinite(state).all():raise FloatingPointError('Nonfinite shared expanded state')
    return state,weights

def select(index,cap,full):
    rows=[r for r in index['cases'] if r['cap']==cap]
    if full:return rows
    lookup={tuple(r['key']):r for r in rows};offsets=list(itertools.product(range(2),repeat=3))
    candidates=[k for k in lookup if all(tuple(np.asarray(k)+d) in lookup for d in offsets)]
    if not candidates:raise ValueError('No multi-axis group')
    # Geometry/count-only selection; never choose based on generated values.
    key=max(candidates,key=lambda k:sum(lookup[tuple(np.asarray(k)+d)]['rows'] for d in offsets))
    return [lookup[tuple(np.asarray(key)+d)] for d in offsets]

def load_observations(rows,chart,expected_binding):
    observations=[];receipts=[]
    for i,r in enumerate(rows):
        folder=INPUT/'tiles'/r['id']/'prepared';receipt=read(folder/'COMPLETE.json')
        if not receipt['passed'] or receipt['tile']!=r['id'] or receipt['binding']!=expected_binding or sha(folder/'arrays.npz')!=receipt['sha256']:
            raise ValueError('Prepared input binding failure: '+r['id'])
        with np.load(folder/'arrays.npz') as f:
            arrays={k:f[k] for k in ('joint','old_joint','wide')}
        if not all(np.isfinite(v).all() for v in arrays.values()):raise ValueError('Invalid observations')
        def t(x,key):return torch.as_tensor(views.scale_channels(x,chart[key]),device='cuda')[None]
        observations.append((t(arrays['old_joint'],'joint'),t(arrays['joint'],'joint'),t(arrays['wide'],'wide')))
        receipts.append({'tile':r['id'],'receipt_sha256':sha(folder/'COMPLETE.json'),'payload_sha256':receipt['sha256']})
        if i%100==0:print(json.dumps({'stage':'conditions','done':i+1,'total':len(rows)}),flush=True)
    return observations,receipts

def addressed_noise(g,cap,seed,stage,draw=0,draws=1):
    rows=[]
    for index in range(draw,draw+draws):
        address=json.dumps(['loa-shared-sixfit-v1',cap,seed,stage,index],separators=(',',':'))
        rng=torch.Generator(device='cuda').manual_seed(int(hashlib.sha256(address.encode()).hexdigest()[:15],16))
        rows.append(torch.randn((1,1,*g.shape),device='cuda',generator=rng))
    value=torch.cat(rows)
    return value if stage=='coarse' else project(value)

def coarse_base(obs):
    return [FieldCondition(old,wide,old.new_zeros((len(old),3)),'wide') for old,new,wide in obs]

def fine_conditions(obs,base,z,cg,expanded):
    result=[]
    for i,(old,new,wide) in enumerate(obs):
        crop=z[cg.crop(i)]
        if expanded:result.append(Condition(new,wide,crop,'sampled'))
        else:result.append(replace(base[i],region='joint',coarse_joint=replicate(crop[:,:,16:32,18:30,18:30]),coarse_wide=crop,coarse_source='sampled'))
    return result

def check_decoding(z,u,cg,fg,chart,expanded):
    cut=(slice(12,36),slice(16,32),slice(16,32)) if expanded else (slice(16,32),slice(18,30),slice(18,30))
    parts=[];mass=0.
    for i in range(len(fg.starts)):
        w=z[cg.crop(i)][(slice(None),slice(None),*cut)].double()
        coarse=torch.exp(w*chart['coarse_logrho']['std'][0]+chart['coarse_logrho']['mean'][0])
        rho=decode_density(coarse,u[fg.crop(i)].double()*chart['fine_residual']['std'][0])
        if not torch.isfinite(rho).all() or not (rho>0).all():raise ValueError('Invalid decoded density')
        mass=max(mass,float(((block_mean(rho)-coarse)/coarse).abs().max()));parts.append(rho)
    comparisons=0;overlap=0.
    for i in range(len(parts)):
        for j in range(i):
            a,b=np.asarray(fg.starts[i]),np.asarray(fg.starts[j]);lo=np.maximum(a,b);hi=np.minimum(a+fg.window,b+fg.window)
            if np.any(hi<=lo):continue
            sl=lambda p:(slice(None),slice(None),*(slice(int(l-q),int(h-q)) for l,h,q in zip(lo,hi,p)))
            overlap=max(overlap,float((parts[i][sl(a)]-parts[j][sl(b)]).abs().max()));comparisons+=1
    if mass>2e-6 or overlap>1e-10 or comparisons!=28:raise ValueError('Mass/overlap gate failed')
    return dict(max_mass_error=mass,overlap_maxabs=overlap,intersection_checks=comparisons)

def source_binding():
    files={Path(__file__).resolve()}
    for module in list(sys.modules.values()):
        filename=getattr(module,'__file__',None)
        if filename and str(filename).startswith(str(ILLUSTRIS)) and str(filename).endswith('.py'):files.add(Path(filename).resolve())
    return {str(p):sha(p) for p in sorted(files)}

@torch.no_grad()
def run(args):
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid') or torch.cuda.device_count()!=1:
        raise RuntimeError('One explicitly allocated GPU required')
    settings.setup(False)
    # Expanded pilot sampling permits PyTorch's native fused attention selection.
    # The older float32 model remains on its supported path; do not benchmark the
    # BF16 transformer under the unrelated old gallery's math-only restriction.
    torch.backends.cuda.enable_flash_sdp(True);torch.backends.cuda.enable_mem_efficient_sdp(True)
    args.out.mkdir(parents=True,exist_ok=False)
    index=read(INPUT/'INDEX.json');index_hash=sha(INPUT/'INDEX.json')
    # The existing prepared-bundle digest convention is deliberately reused.
    import loa_survey_runner_20261016 as runner
    expected=runner.digest(index);sources=source_binding()
    binding=dict(index_sha256=index_hash,source=sources,checkpoints=index['checkpoints'],mode=args.mode,
                 cap=args.cap,seed=args.seed,job=os.environ['SLURM_JOB_ID'],truth_opened=False,
                 retraining=False,posterior_qualified=False,tf32=False,draws=args.draws,expanded_autocast='bfloat16',fused_attention_allowed=True,
                 gpu_name=torch.cuda.get_device_name(0),gpu_total_gib=torch.cuda.get_device_properties(0).total_memory/2**30,
                 nfe=128 if args.mode=='smoke' else 2)
    write(args.out/'BINDING.json',binding);results=[];started=time.monotonic()
    def progress(label):
        def report(step,state):
            elapsed=time.monotonic()-started
            if elapsed>args.seconds:raise TimeoutError('Engineering time budget exhausted; partial receipts retained')
            if step%16==0 or binding['nfe']==2:print(json.dumps({'stage':label,'heun_step':step,'elapsed_seconds':elapsed}),flush=True)
        return report
    seeds=(17,29) if args.seed is None else (args.seed,)
    caps=('NGC','SGC') if args.cap is None else (args.cap,)
    for seed in seeds:
        models,hashes,normalizer=evaluate.load_models(seed,13312,26624)
        for stage,digest in hashes.items():
            if digest!=index['checkpoints'][f'{"fine" if stage=="fine" else "coarse"}_{seed}']['sha256']:raise ValueError('Checkpoint drift')
        chart=views.load_chart(normalizer)
        for cap in caps:
            rows=select(index,cap,args.mode=='benchmark');mid=[r['midpoint'] for r in rows];cg=Geometry.make(mid,'coarse')
            tick=time.monotonic();obs,inputs=load_observations(rows,chart,expected)
            obs=[tuple(v.expand(args.draws,*v.shape[1:]) for v in item) for item in obs]
            base=coarse_base(obs)
            write(args.out/f'INPUTS_{cap}_{seed}.json',inputs)
            preparation=time.monotonic()-tick;torch.cuda.reset_peak_memory_stats();tick=time.monotonic()
            nfe=binding['nfe'];z,cw=integrate(models['coarse'],'coarse',base,addressed_noise(cg,cap,seed,'coarse',draws=args.draws),cg,nfe,progress(f'{cap}_{seed}_coarse'))
            torch.cuda.synchronize();coarse_seconds=time.monotonic()-tick
            coarse_peak=torch.cuda.max_memory_allocated()/2**30
            coarse_sha=hashlib.sha256(z.cpu().numpy().tobytes()).hexdigest()
            for family in FAMILIES:
                expanded=family!='legacy';fg=geometry(mid,expanded);fine=fine_conditions(obs,base,z,cg,expanded)
                model=models['fine'];model_id=f'legacy_v_{seed}';parity=None
                if expanded:
                    model_id=f'{family}_{"x" if family=="unet" else "v"}_{seed}';ck=index['checkpoints'][model_id]
                    if sha(ck['path'])!=ck['sha256']:raise ValueError('Expanded checkpoint drift')
                    state=torch.load(ck['path'],map_location='cpu',weights_only=False)
                    if state['step']!=9984 or state['binding']['sources']!=training_sources() or state['binding']['normalizer']!=normalizer:raise ValueError('Expanded model binding mismatch')
                    net=FineModel(family,'x' if family=='unet' else 'v').cuda().eval();net.load_state_dict(state['ema']);del state
                    model=ExpandedAdapter(net)
                    if args.mode in ('smoke','benchmark'):
                        single=geometry([mid[0]],True);noise=addressed_noise(single,cap,seed,'expanded-parity',draws=args.draws)
                        direct=sample(net,fine[0],noise,4)
                        coupled,_=integrate_expanded(net,[fine[0]],noise,single,4)
                        parity=float((direct-coupled).abs().max());assert parity<1e-6
                        del direct,coupled,noise
                torch.cuda.reset_peak_memory_stats();tick=time.monotonic()
                initial=addressed_noise(fg,cap,seed,'expanded-fine' if expanded else 'legacy-fine',draws=args.draws)
                if expanded:u,fw=integrate_expanded(net,fine,initial,fg,nfe,progress(cap+'_'+model_id))
                else:u,fw=integrate(model,'fine',fine,initial,fg,nfe,progress(cap+'_'+model_id))
                torch.cuda.synchronize();elapsed=time.monotonic()-tick
                item=dict(cap=cap,seed=seed,model=model_id,windows=len(rows),draws=args.draws,nfe=nfe,
                          preparation_seconds=preparation,coarse_seconds=coarse_seconds,fine_seconds=elapsed,
                          coarse_latent_sha256=coarse_sha,single_window_parity_maxabs=parity,
                          peak_gpu_gib=torch.cuda.max_memory_allocated()/2**30,coarse_peak_gpu_gib=coarse_peak,
                          fine_shape=list(map(int,fg.shape)),coarse_shape=list(map(int,cg.shape)),
                          posterior_output=False)
                if args.mode=='smoke':item.update(check_decoding(z,u,cg,fg,chart,expanded))
                results.append(item);write(args.out/(model_id+'_'+cap+'.json'),item);print(json.dumps(item),flush=True)
                del u,fw,fine,model,initial
                if expanded:del net
                gc.collect();torch.cuda.empty_cache()
            del obs,base,z,cw;gc.collect();torch.cuda.empty_cache()
        del models;gc.collect();torch.cuda.empty_cache()
    if source_binding()!=sources or sha(INPUT/'INDEX.json')!=index_hash:raise ValueError('Source/input changed during test')
    write(args.out/'COMPLETE.json',dict(passed=True,results=results,seconds=time.monotonic()-started,binding_sha256=sha(args.out/'BINDING.json'),posterior_qualified=False))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['smoke','benchmark'],required=True);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--cap',choices=['NGC','SGC']);p.add_argument('--seed',type=int,choices=[17,29]);p.add_argument('--seconds',type=float,default=3000)
    p.add_argument('--draws',type=int,choices=[1,4],default=1);a=p.parse_args();run(a)
