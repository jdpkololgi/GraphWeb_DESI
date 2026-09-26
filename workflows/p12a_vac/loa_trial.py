"""Versioned Loa halo48 canary: observer-only input preparation and provisional VAC."""
import argparse,sys,os,json,hashlib,importlib.util,dataclasses,time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import fitsio
import healpy as hp
import torch
ILL=Path(os.environ.get('ILLUSTRIS_ROOT','/global/u2/d/dkololgi/TNG/Illustris'))
OBS=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ILL))
spec=importlib.util.spec_from_file_location('observation_patch',OBS/'workflows/catalog/p12a_observation_patch.py');obs=importlib.util.module_from_spec(spec);spec.loader.exec_module(obs)
from workflows.abacus_tweb.p3br_build_random_response import add_random_file,normalized_map,angular_boundary_distance,galactic_cap,photsys_code,NSIDE
from workflows.abacus_tweb.p3a_build_canonical_fields import grid_from_xyz
from workflows.abacus_tweb.p6_field_patch_utils import FieldPatch
from workflows.abacus_tweb import p8_train_unet_patch as u
from workflows.abacus_tweb.p8_deterministic_common import increments_to_eigenvalues,unscale_increments
from workflows.sbi.p12a_blind_inference import reconstruct_fmpe
from workflows.sbi.p12_train_base_response_fmpe import sample_posterior,theta_to_eigenvalues
from workflows.sbi.p12_export_unet_summaries import ntilde_at_rows
B=Path('/pscratch/sd/d/dkololgi/abacus/p10_multiphase');C=B/'p12a_halo48_candidate_20260924_v1'
def predict(model,patch,ck):
 x,p=u.model_inputs(patch,ck['normalization'],'cuda')
 with torch.inference_mode():scaled=model.head(model.sample_latent(x,p)).cpu().numpy()
 return increments_to_eigenvalues(unscale_increments(scaled,ck['scaler'])).astype('f4')

def read(p):return json.loads(Path(p).read_text())
def digest(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,d):
 t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n');t.replace(p)
def response_at(xyz,grid,angular,angles,selection):
 # P3b stores distances at voxel centres and samples with floor((x-origin)/cell).
 index=np.floor((xyz-np.asarray(grid['origin_mpc']))/grid['cell_mpc']).astype('i8')
 xyz=np.asarray(grid['origin_mpc'])+(index+.5)*grid['cell_mpc'];radius=np.linalg.norm(xyz,axis=1)
 pix=hp.vec2pix(NSIDE,*xyz.T,nest=False)
 z=np.interp(radius,selection['cosmology']['radius_grid_mpc'],selection['cosmology']['redshift_grid'])
 support=angular[pix]&(z>=.1)&(z<.6)&~((z>=.585)&(z<.595))
 edges=np.interp([.1,.6,.585,.595],selection['cosmology']['redshift_grid'],selection['cosmology']['radius_grid_mpc'])
 distance=np.minimum(angles[pix].astype('f8')*radius,np.min(abs(radius[:,None]-edges),axis=1)).astype('f4');distance*=support
 return distance,support

def contiguous_random_file(counts,record):
 accepted=rejected=maskbits=0
 with fitsio.FITS(record['path']) as f:
  n=f[1].get_nrows()
  for start in range(0,n,250000):
   # Read complete contiguous records; project columns only in memory.
   block=f[1][start:min(start+250000,n)]
   good=np.asarray(block['GOODHARDLOC'],dtype=bool);rejected+=int((~good).sum());maskbits+=int(np.count_nonzero(block['MASKBITS']))
   if not good.any():continue
   ra=np.asarray(block['RA'][good],dtype='f8');dec=np.asarray(block['DEC'][good],dtype='f8')
   pix=hp.ang2pix(NSIDE,ra,dec,lonlat=True,nest=False);domain=galactic_cap(ra,dec).astype('i8')*2+photsys_code(block['PHOTSYS'][good]).astype('i8')
   counts+=np.bincount(domain*counts.shape[1]+pix,minlength=counts.size).reshape(counts.shape);accepted+=int(good.sum())
 return dict(rows=n,accepted_rows=accepted,goodhardloc_rejected=rejected,maskbits_nonzero=maskbits,maskbits_filter_applied=False,read_mode='contiguous whole-row blocks; same registered selection and counts')

def prepare(root):
 if (root/'INPUTS_READY.json').exists():raise FileExistsError('inputs already complete')
 root.mkdir(parents=True,exist_ok=True);(root/'random_maps').mkdir(exist_ok=True)
 registry=read(ILL/'configs/p10_response_sources_v1.json')['desi_candidate'];source=registry['data']['full'];path=Path(source['path'])
 if digest(path)!=source['sha256']:raise ValueError('Loa full catalogue hash changed')
 if (root/'catalogue.npz').exists() and (root/'galaxy_angular_support.npy').exists():
  cached=np.load(root/'catalogue.npz');xyz=cached['xyz'];cap=cached['cap'];shell=cached['shell'];d=cached['targetid'];quality_rows=None
  if len(np.unique(d))!=len(d):raise ValueError('duplicate cached TARGETID')
 else:
  cols=['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE'];parts=[];angular_counts=np.zeros(hp.nside2npix(256),dtype='i8');quality_rows=0
  with fitsio.FITS(path) as f:
   for start in range(0,f[1].get_nrows(),500000):
    d=f[1].read(rows=range(start,min(start+500000,f[1].get_nrows())),columns=cols);good=obs.successful_rows(d);d=d[good];quality_rows+=len(d)
    angular_counts+=np.bincount(hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True),minlength=len(angular_counts))
    z=d['Z_not4clus'];d=d[(z>=.1)&(z<.6)&~((z>=.585)&(z<.595))];parts.append(d)
  d=np.concatenate(parts);del parts
  if len(np.unique(d['TARGETID']))!=len(d):raise ValueError('duplicate TARGETID')
  xyz=obs.observer_xyz(d['RA'],d['DEC'],d['Z_not4clus']);cap=galactic_cap(d['RA'],d['DEC']);shell=np.searchsorted([.15,.25,.35,.45,.55],d['Z_not4clus'],side='right')-1;shell[(shell<0)|(shell>3)]=-1
  np.savez(root/'catalogue.npz',targetid=d['TARGETID'],ra=d['RA'],dec=d['DEC'],z=d['Z_not4clus'],xyz=xyz,cap=cap,shell=shell)
  np.save(root/'galaxy_angular_support.npy',angular_counts>0)
 grids={name:grid_from_xyz(xyz[cap==i],5.,40.).as_dict() for i,name in [(0,'SGC'),(1,'NGC')]}
 def one(pair):
  i,record=pair;p=root/f'random_maps/{i:02d}.npz';q=root/f'random_maps/{i:02d}.json'
  if p.exists() and q.exists():
   receipt=read(q)
   if receipt['source_sha256']!=record['sha256'] or receipt['map_sha256']!=digest(p):raise ValueError('random map binding changed')
   return np.load(p)['counts']
  before=Path(record['path']).stat();counts=np.zeros((4,hp.nside2npix(NSIDE)),dtype='i8');audit=contiguous_random_file(counts,record);after=Path(record['path']).stat()
  if (before.st_size,before.st_mtime_ns)!=(after.st_size,after.st_mtime_ns) or before.st_size!=record['bytes']:raise ValueError('random source changed')
  np.savez_compressed(p,counts=counts);save(q,dict(source_sha256=record['sha256'],map_sha256=digest(p),audit=audit));print('random',i,flush=True);return counts
 with ThreadPoolExecutor(max_workers=4) as pool:counts=sum(pool.map(one,enumerate(registry['full_random'])))
 angular=normalized_map(counts);np.savez_compressed(root/'random_angular.npz',support=angular['support'],domain=angular['domain'],response=angular['angular_response'])
 angles={i:angular_boundary_distance(angular['support']&((angular['domain']//2)==i)) for i in [0,1]};np.savez(root/'boundary_angles.npz',SGC=angles[0],NGC=angles[1])
 # Truth-free preselection: median-occupancy and strongest-boundary core per cap/shell.
 cases=[]
 selection=read(B/'training_contract/transforms/field/selection_manifest.json')
 for i,name in [(0,'SGC'),(1,'NGC')]:
  grid=grids[name];rows=np.flatnonzero((cap==i)&(shell>=0));core=np.floor((xyz[rows]-grid['origin_mpc'])/5/20).astype('i8');unique,inverse,n=np.unique(core,axis=0,return_inverse=True,return_counts=True)
  centres=np.asarray(grid['origin_mpc'])+(unique*20+10)*5;zc=np.interp(np.linalg.norm(centres,axis=1),selection['cosmology']['radius_grid_mpc'],selection['cosmology']['redshift_grid']);ss=np.searchsorted([.15,.25,.35,.45,.55],zc,side='right')-1
  pix=hp.vec2pix(NSIDE,*centres.T);edge=angles[i][pix]*np.linalg.norm(centres,axis=1)
  for sh in range(4):
   eligible=np.flatnonzero((ss==sh)&(n>=16))
   if len(eligible)<2:raise ValueError('insufficient canary strata')
   med=eligible[np.argsort(n[eligible],kind='stable')[len(eligible)//2]];remaining=eligible[eligible!=med];boundary=remaining[np.lexsort((remaining,edge[remaining]))[0]]
   for kind,j in [('typical',med),('boundary',boundary)]:cases.append(dict(cap=i,shell=sh,kind=kind,core_start=(unique[j]*20).tolist(),core_stop=np.minimum(unique[j]*20+20,grid['shape']).tolist(),active_rows=int(n[j])))
 save(root/'INPUTS_READY.json',dict(schema='loa-p12a-canary-inputs-v1',source=source,quality_rows=quality_rows,context_rows=len(d),active_rows=int((shell>=0).sum()),grids=grids,cases=cases,selection_path=str(B/'training_contract/transforms/field/selection_manifest.json'),source_registry_sha256=digest(ILL/'configs/p10_response_sources_v1.json'),arrays_sha256={p.name:digest(p) for p in [root/'catalogue.npz',root/'galaxy_angular_support.npy',root/'random_angular.npz',root/'boundary_angles.npz']},provisional=True))
 print('INPUTS READY',len(d),len(cases),flush=True)

def infer(root):
    parity=read(root/'RESPONSE_PARITY.json')
    if not parity['distance_exact'] or not parity['support_exact']:raise ValueError('response parity gate failed')
    ready=read(root/'INPUTS_READY.json')
    for name,sha in ready['arrays_sha256'].items():
        if digest(root/name)!=sha:raise ValueError('input hash mismatch: '+name)
    if (root/'VAC_COMPLETE.json').exists():raise FileExistsError('VAC already complete')
    selection=read(ready['selection_path']);p3=read(B/'ph006/p3_fields/field_manifest.json');schema=read(p3['frozen_schema']);spline=read(p3['ntilde_spline'])
    for grid in ready['grids'].values():
        if grid['padding_mpc']!=schema['grid']['padding_mpc']:raise ValueError('grid padding differs from frozen schema')
    catalogue=np.load(root/'catalogue.npz');xyz=catalogue['xyz'];cap=catalogue['cap'];shell=catalogue['shell'];z=catalogue['z']
    angular=np.load(root/'galaxy_angular_support.npy');random=np.load(root/'random_angular.npz');angles=np.load(root/'boundary_angles.npz')
    sr=read(C/'summaries/ph006/OOF_SUMMARY_COMPLETE.json');ck=torch.load(sr['checkpoint'],map_location='cuda',weights_only=False)
    if digest(sr['checkpoint'])!=sr['checkpoint_sha256']:raise ValueError('encoder changed')
    torch.set_num_threads(8);torch.backends.cudnn.allow_tf32=True;torch.backends.cuda.matmul.allow_tf32=False
    if digest(C/'posterior/fmpe_estimator.pt')!=read(C/'posterior/P12A_COMPLETE.json')['checkpoint_sha256']:raise ValueError('posterior checkpoint changed')
    model=u.UPatch().cuda().eval();model.load_state_dict(ck['state_dict']);posterior,pck=reconstruct_fmpe(C/'posterior/fmpe_estimator.pt','cuda')
    train=read(C/'dataset/P12A_DATASET_READY.json')
    with np.load(train['training']['path']) as f:lo=f['context'].min(axis=0);hi=f['context'].max(axis=0)
    outputs=[];draw_files=[];qa=[];seen=set()
    (root/'shards').mkdir(exist_ok=True)
    for k,case in enumerate(ready['cases']):
        name='NGC' if case['cap'] else 'SGC';grid=ready['grids'][name];origin=np.asarray(grid['origin_mpc']);start=np.asarray(case['core_start']);stop=np.asarray(case['core_stop']);shape=np.asarray(grid['shape'])
        frac=(xyz-origin)/5-.5;voxel=np.floor((xyz-origin)/5).astype('i8')
        rows=np.flatnonzero((cap==case['cap'])&(shell>=0)&np.all((voxel>=start)&(voxel<stop),axis=1))
        ids=catalogue['targetid'][rows]
        if seen.intersection(ids.tolist()):raise ValueError('output core ownership collision')
        seen.update(ids.tolist());cs=np.maximum(((start-48)//8)*8,0);ce=np.minimum(((stop+48+7)//8)*8,shape)
        context=np.flatnonzero((cap==case['cap'])&np.all((frac>=cs-1)&(frac<ce),axis=1))
        fields=obs.rebuild_fields(xyz[context],shell[context],grid,cs,ce,angular,256,schema,spline,selection,0,name)
        if not all(np.isfinite(fields[n]).all() for n in u.CHANNELS):raise ValueError('nonfinite input fields')
        patch=FieldPatch(k,-1,case['cap'],tuple(u.CHANNELS),np.stack([fields[n] for n in u.CHANNELS]),cs,ce,start,stop,tuple(slice(int(a),int(b)) for a,b in zip(start-cs,stop-cs)),rows,frac[rows],frac[rows]-cs,start-cs,ce-stop)
        base=predict(model,patch,ck)
        distance,support=response_at(xyz[rows],grid,random['support']&((random['domain']//2)==case['cap']),angles[name],selection)
        nt=ntilde_at_rows(selection,cap[rows],z[rows].astype('f4'))
        x=np.column_stack([base,z[rows].astype('f4'),np.log(nt),cap[rows],np.log1p(distance)]).astype('f4')
        outside=np.any((x<lo)|(x>hi),axis=1);eligible=support&np.isfinite(x).all(axis=1)
        draws=np.full((len(rows),512,3),np.nan,dtype='f4')
        torch.manual_seed(20260925+k);torch.cuda.manual_seed_all(20260925+k)
        if eligible.any():
            scaled=((x[eligible]-pck['context_mean'])/pck['context_std']).astype('f4');sample=sample_posterior(posterior,scaled,512,128,'cuda')
            draws[eligible]=theta_to_eigenvalues(sample*np.asarray(pck['theta_std'])+np.asarray(pck['theta_mean'])).astype('f4')
            if not np.isfinite(draws[eligible]).all() or np.any(np.diff(draws[eligible],axis=-1)<0):raise ValueError('invalid posterior draws')
        dtype=[('TARGETID','i8'),('RA','f8'),('DEC','f8'),('Z','f8'),('CAP','i1'),('CORE_ID','i4'),('QUALITY','u2'),('SUPPORTED','?'),('BASE_EIGENVALUES','f4',(3,)),('EIGENVALUE_MEAN','f4',(3,)),('EIGENVALUE_Q05','f4',(3,)),('EIGENVALUE_Q16','f4',(3,)),('EIGENVALUE_Q50','f4',(3,)),('EIGENVALUE_Q84','f4',(3,)),('EIGENVALUE_Q95','f4',(3,)),('P_WEB','f4',(4,)),('BOUNDARY_MPC','f4'),('NTILDE_MPC3','f4')]
        out=np.zeros(len(rows),dtype=dtype)
        for n,v in [('TARGETID',ids),('RA',catalogue['ra'][rows]),('DEC',catalogue['dec'][rows]),('Z',z[rows]),('CAP',cap[rows]),('CORE_ID',k),('SUPPORTED',eligible),('BASE_EIGENVALUES',base),('BOUNDARY_MPC',distance),('NTILDE_MPC3',nt)]:out[n]=v
        out['QUALITY']=np.uint16(32)+(~eligible).astype('u2')+2*(shell[rows]==3).astype('u2')+4*(distance<10.345846881466155).astype('u2')+8*(distance<20.69169376293231).astype('u2')+16*outside.astype('u2')
        for n in ['EIGENVALUE_MEAN','EIGENVALUE_Q05','EIGENVALUE_Q16','EIGENVALUE_Q50','EIGENVALUE_Q84','EIGENVALUE_Q95','P_WEB']:out[n]=np.nan
        if eligible.any():
            out['EIGENVALUE_MEAN'][eligible]=draws[eligible].mean(axis=1);q=np.quantile(draws[eligible],[.05,.16,.5,.84,.95],axis=1)
            for j,n in enumerate(['EIGENVALUE_Q05','EIGENVALUE_Q16','EIGENVALUE_Q50','EIGENVALUE_Q84','EIGENVALUE_Q95']):out[n][eligible]=q[j]
            classes=np.sum(draws[eligible]>.2,axis=-1)
            out['P_WEB'][eligible]=np.stack([(classes==j).mean(axis=1) for j in range(4)],axis=1)
            if not np.allclose(out['P_WEB'][eligible].sum(axis=1),1):raise ValueError('class probabilities not normalized')
        file=root/f'shards/core_{k:02d}_draws.npz';np.savez_compressed(file,TARGETID=ids,eigenvalue_draws=draws);draw_files.append(dict(path=str(file),sha256=digest(file)))
        outputs.append(out);qa.append(dict(**case,rows=len(rows),supported=int(eligible.sum()),outside_training_envelope=int(outside.sum()),context_rows=len(context),raw_input_count_sum=float(fields['counts'].sum(dtype='f8')),median_base= np.median(base,axis=0).tolist()))
        print('core',k,'rows',len(rows),'supported',eligible.sum(),'outside',outside.sum(),flush=True)
    out=np.concatenate(outputs);out.sort(order='TARGETID');file=root/'DESI_LOA_P12A_HALO48_CANARY_VAC.fits'
    fitsio.write(file,out,extname='ENV_POSTERIOR',header={'PROVIS':True,'MODEL':'P12A_H48','NSAMPLE':512,'LTHRESH':.2,'RSMOOTH':7.,'RUNIT':'Mpc/h','ZTARGET':.2},clobber=False)
    reread=fitsio.read(file);assert np.array_equal(reread['TARGETID'],out['TARGETID'])
    save(root/'VAC_COMPLETE.json',dict(schema='desi-loa-p12a-halo48-canary-vac-v1',technical_complete=True,science_release_ready=False,scope='16 geometry-selected diagnostic cores, not full footprint or random galaxy sample',rows=len(out),supported_rows=int(out['SUPPORTED'].sum()),catalogue=str(file),catalogue_sha256=digest(file),draw_shards=draw_files,cases=qa,inputs_marker_sha256=digest(root/'INPUTS_READY.json'),posterior_checkpoint_sha256=digest(C/'posterior/fmpe_estimator.pt'),quality_bits={'1':'unsupported; posterior null','2':'sparse shell z>=0.45','4':'boundary distance <R/h','8':'boundary distance <2R/h','16':'one or more posterior features outside training min/max envelope','32':'provisional diagnostic only'},class_order=['void','sheet','filament','knot'],job=os.environ['SLURM_JOB_ID']))
    print('VAC COMPLETE',len(out),flush=True)


def finalize_geometry(root):
 ready=read(root/'INPUTS_READY.json');catalogue=np.load(root/'catalogue.npz');xyz=catalogue['xyz'];cap=catalogue['cap'];shell=catalogue['shell']
 padding=read(B/'ph006/p3_fields/p3_field_schema_v1.json')['grid']['padding_mpc']
 grids={name:grid_from_xyz(xyz[cap==i],5.,padding).as_dict() for i,name in [(0,'SGC'),(1,'NGC')]}
 cached=np.load(root/'boundary_angles.npz');angles={0:cached['SGC'],1:cached['NGC']}
 # Truth-free preselection: median-occupancy and strongest-boundary core per cap/shell.
 cases=[]
 selection=read(B/'training_contract/transforms/field/selection_manifest.json')
 for i,name in [(0,'SGC'),(1,'NGC')]:
  grid=grids[name];rows=np.flatnonzero((cap==i)&(shell>=0));core=np.floor((xyz[rows]-grid['origin_mpc'])/5/20).astype('i8');unique,inverse,n=np.unique(core,axis=0,return_inverse=True,return_counts=True)
  centres=np.asarray(grid['origin_mpc'])+(unique*20+10)*5;zc=np.interp(np.linalg.norm(centres,axis=1),selection['cosmology']['radius_grid_mpc'],selection['cosmology']['redshift_grid']);ss=np.searchsorted([.15,.25,.35,.45,.55],zc,side='right')-1
  pix=hp.vec2pix(NSIDE,*centres.T);edge=angles[i][pix]*np.linalg.norm(centres,axis=1)
  for sh in range(4):
   eligible=np.flatnonzero((ss==sh)&(n>=16))
   if len(eligible)<2:raise ValueError('insufficient canary strata')
   med=eligible[np.argsort(n[eligible],kind='stable')[len(eligible)//2]];remaining=eligible[eligible!=med];boundary=remaining[np.lexsort((remaining,edge[remaining]))[0]]
   for kind,j in [('typical',med),('boundary',boundary)]:cases.append(dict(cap=i,shell=sh,kind=kind,core_start=(unique[j]*20).tolist(),core_stop=np.minimum(unique[j]*20+20,grid['shape']).tolist(),active_rows=int(n[j])))
 ready['grids']=grids;ready['cases']=cases;ready['geometry_source_sha256']=digest(Path(__file__))
 previous=root/'INPUTS_READY_INITIAL_GEOMETRY.json'
 if previous.exists():raise FileExistsError(previous)
 (root/'INPUTS_READY.json').rename(previous);save(root/'INPUTS_READY.json',ready)
 print('FINAL GEOMETRY READY',len(cases),padding,flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','finalize','infer']);p.add_argument('--root',type=Path,required=True);a=p.parse_args()
 if a.stage=='prepare':prepare(a.root)
 elif a.stage=='finalize':finalize_geometry(a.root)
 else:infer(a.root)
