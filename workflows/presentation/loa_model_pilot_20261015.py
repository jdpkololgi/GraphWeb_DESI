"""Targetless real-Loa pilot: CPU preparation, then frozen field draws on one GPU.

No training, simulation truth, automatic allocation, or full-cap claim. Existing
products are read-only. The new products share owned cores, observation recipes,
coarse draws within each seed, and a declared finite local/wide tidal closure.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '4'
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import argparse
import hashlib
import json
from pathlib import Path
import socket
import sys
import time

ILLUSTRIS = Path('/global/u2/d/dkololgi/TNG/Illustris')
sys.path.insert(0, str(ILLUSTRIS))
ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_model_pilot_20261015_v1')
ATLAS = Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1')
VAC = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')
PHOT = Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def save_json(path, obj):
    with Path(path).open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)


def require_compute():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Payload work requires an allocated compute node')


def prepare():
    require_compute()
    import fitsio
    import h5py
    import numpy as np
    from workflows.sbi import e2e_loa_field_conditions as lc
    from workflows.sbi import e2e_wholefield_pilot_data as wd
    from workflows.sbi import e2e_coupled_contract as c
    from workflows.sbi.e2e_loa_field_atlas import tile_keys
    start = time.monotonic()
    ROOT.mkdir(parents=True, exist_ok=False)
    import torch
    from workflows.sbi.e2e_wholefield_pilot_train import sources
    checkpoints = {}
    for architecture, prediction in [('unet','x'),('transformer','v')]:
        for seed in (17,29):
            modelid = f'{architecture}_{prediction}_{seed}'
            path = wd.ROOT/'fits'/modelid/'CHECKPOINT_009984.pt'
            state = torch.load(path,map_location='cpu',weights_only=False)
            if state['binding']['sources'] != sources() or state['step'] != 9984:
                raise ValueError('Frozen model source mismatch before inference')
            checkpoints[modelid] = sha(path)
            del state
    plan = read(ATLAS/'PLAN.json')
    # Selection depends only on observed row counts, not properties or outcomes.
    cases = [max((r for r in plan['tiles'] if r['cap'] == cap), key=lambda r: (r['rows'], r['id']))
             for cap in ('NGC', 'SGC')]
    source = {'runner': sha(__file__), 'atlas_plan': sha(ATLAS/'PLAN.json'), 'expanded_checkpoints':checkpoints}
    for cap in ('NGC', 'SGC'):
        receipt = read(ATLAS/f'{cap}_CACHE.json')
        record = receipt['outputs'][0]
        if sha(record['path']) != record['sha256']:
            raise ValueError('Observation cache hash mismatch')
        source[cap+'_cache_receipt'] = sha(ATLAS/f'{cap}_CACHE.json')
    if sha(ATLAS/'vac_properties.npz') != plan['properties_sha256']:
        raise ValueError('Saved property catalogue changed')
    with np.load(ATLAS/'vac_properties.npz') as f:
        xyz, capids, ids = f['xyz'], f['CAP'], f['TARGETID']
    selected = []
    for case in cases:
        cap = case['cap']; grid = plan['grids'][cap]; mid = np.asarray(case['midpoint'])
        rows = np.flatnonzero(capids == case['cap_id'])
        keys = tile_keys(xyz[rows], grid)
        rows = rows[(keys == np.asarray(case['key'])).all(1)]
        assert len(rows) == case['rows']
        case['grid'] = grid
        case['core_origin_mpc_h'] = (np.asarray(grid['origin_mpc_h'])+(mid-[32, 16, 16])*3.383).tolist()
        position = (xyz[rows]-np.asarray(case['core_origin_mpc_h']))/6.766-.5
        inside = ((position >= 0) & (position <= [31, 15, 15])).all(1)
        case['excluded_interpolation_edge'] = int((~inside).sum())
        rows = rows[inside]; selected.append(rows)
        response = lc.load_response(cap)
        with h5py.File(ATLAS/f'{cap}_observation_cache.h5', 'r') as cache:
            old = lc.build(grid, mid, cap, cache['counts_raw'], cache['wide'], response)
            expanded = wd.build_fine_observations(grid, mid, cap, cache['counts_raw'], response)
        error = float(np.max(abs(expanded[(slice(None), *wd.OLD_CROP)]-old['joint'])))
        if error > 2e-5:
            raise ValueError('Expanded/legacy central observation mismatch')
        directory = ROOT/cap; directory.mkdir()
        np.savez_compressed(directory/'conditions.npz', joint=expanded, old_joint=old['joint'], wide=old['wide'])
        np.savez_compressed(directory/'positions.npz', xyz=xyz[rows], targetid=ids[rows],
                            position=(xyz[rows]-np.asarray(case['core_origin_mpc_h']))/6.766-.5)
        case.update(central_channel_maxabs=error, retained_rows=len(rows),
                    condition_sha256=sha(directory/'conditions.npz'))
        print(json.dumps({'stage': 'conditions', 'cap': cap, 'rows': len(rows), 'seconds': time.monotonic()-start}), flush=True)
    chosen = np.concatenate(selected)
    cols = ['TARGETID','RA','DEC','Z','SUPPORTED','CIGALE_MATCHED','CIGALE_AMBIGUOUS','CIGALE_SPECTYPE',
            'MASS_CG_15','SFR_CG_15','AGE_CG_15','AV_CG_15','FAGN_CG_15','P_WEB']
    with fitsio.FITS(VAC) as f:
        a = f[1].read(rows=np.sort(chosen), columns=cols)
    byid = np.argsort(a['TARGETID']); sortedids = a['TARGETID'][byid]
    assert (np.diff(sortedids) > 0).all()
    r = np.full(len(a), np.nan); matched = np.zeros(len(a), dtype='i4')
    with fitsio.FITS(PHOT) as f:
        for startrow in range(0, f[1].get_nrows(), 200000):
            b = f[1].read(rows=np.arange(startrow, min(startrow+200000, f[1].get_nrows())),
                          columns=['TARGETID','RA','DEC','Z_not4clus','FLUX_R','MW_TRANSMISSION_R','PHOTSYS'])
            at = np.searchsorted(sortedids, b['TARGETID'])
            good = (at < len(a)) & (sortedids[np.minimum(at, len(a)-1)] == b['TARGETID'])
            ii = byid[at[good]]; b = b[good]
            for x,y in [('RA','RA'),('DEC','DEC'),('Z_not4clus','Z')]:
                np.testing.assert_array_equal(a[y][ii], b[x])
            flux = b['FLUX_R']/b['MW_TRANSMISSION_R']
            ps = np.char.strip(b['PHOTSYS'].astype('U'))
            valid = np.isfinite(flux) & (flux > 0) & np.isin(ps, ['N','S'])
            r[ii[valid]] = 22.5-2.5*np.log10(flux[valid])-np.where(ps[valid] == 'N', .04, 0)
            matched[ii] += 1
    assert np.all(matched == 1)
    for case in cases:
        directory = ROOT/case['cap']
        with np.load(directory/'positions.npz') as f:
            targetid = f['targetid']
        at = byid[np.searchsorted(sortedids, targetid)]
        np.testing.assert_array_equal(a['TARGETID'][at], targetid)
        np.savez_compressed(directory/'properties.npz', r_bright=r[at], **{name:a[name][at] for name in cols})
        case['z_range'] = [float(a['Z'][at].min()), float(a['Z'][at].max())]
        case['properties_sha256'] = sha(directory/'properties.npz')
    # Large observed parent counts establish coverage, not an extrapolation claim.
    save_json(ROOT/'PREPARED.json', {'passed':True, 'cases':cases, 'sources':source,
        'parent_counts':plan['counts'], 'supported_redshift':plan['science_redshift'],
        'selection':'highest observed row-count owned tile per cap; throughput/property-rich pilot, not representative survey sample',
        'truth_opened':False, 'seconds':time.monotonic()-start, 'job':os.environ['SLURM_JOB_ID'],
        'host':socket.gethostname(), 'vac':str(VAC), 'photometry':str(PHOT)})
    print('CPU_PREPARATION_COMPLETE', flush=True)


def infer():
    require_compute()
    import numpy as np
    import torch
    from dataclasses import replace
    from scipy.ndimage import map_coordinates
    from workflows.sbi import e2e_coupled_contract as c
    from workflows.sbi import e2e_coupled_views as views
    from workflows.sbi import e2e_coupled_operators as op
    from workflows.sbi import e2e_footprint_tidal_operator as tidal
    from workflows.sbi.e2e_cfm_pilot_evaluate import load_models
    from workflows.sbi.e2e_coupled_benchmark_models import FieldCondition, sample as oldsample
    from workflows.sbi.e2e_vdm_context_models import replicate, decode_density
    from workflows.sbi.e2e_wholefield_pilot_models import FineModel, Condition, sample
    from workflows.sbi.e2e_wholefield_pilot_train import setup, sources
    from workflows.sbi.e2e_wholefield_pilot_data import ROOT as FITROOT, SHAPE, OLD_CROP, COARSE_CROP
    setup(17)
    prep = read(ROOT/'PREPARED.json')
    if prep['sources']['runner'] != sha(__file__):
        raise ValueError('Runner changed after preparation')
    start = time.monotonic(); all_records = []
    cs = (slice(32,64),slice(24,40),slice(24,40))
    oldcs = (slice(16,48),slice(16,32),slice(16,32))
    # Exact free-space core operator, unit density beyond the common wide cube.
    kernels = {}
    for key, shape, offset in [('wide',(192,)*3,(80,88,88)),('new',SHAPE,(32,24,24)),('old',(64,48,48),(16,16,16))]:
        k = tidal.kernel(shape, (32,16,16), offset=offset)
        kernels[key] = (k.shape[1:], [torch.fft.rfftn(torch.as_tensor(v,device='cuda')) for v in k])
    @torch.no_grad()
    def tensor_core(rho, wide, kind):
        lifted = op.lift(wide)
        cut = (slice(48,144),slice(64,128),slice(64,128)) if kind == 'new' else (slice(64,128),slice(72,120),slice(72,120))
        residual = rho-lifted[cut]
        result = np.zeros((32,16,16,6))
        for value, key in [(lifted-1,'wide'),(residual,kind)]:
            shape, ks = kernels[key]
            spectrum = torch.fft.rfftn(torch.as_tensor(value,device='cuda'),s=shape)
            for i,k in enumerate(ks):
                result[...,i] += torch.fft.irfftn(spectrum*k,s=shape)[:32,:16,:16].cpu().numpy()
        rc = rho[cs if kind == 'new' else oldcs]
        error = float(np.max(abs(result[...,[0,3,5]].sum(-1)-(rc-1))))
        if not np.isfinite(result).all() or error > 2e-8*max(1.,float(abs(rc-1).max())):
            raise ValueError('Core free-space trace check failed')
        return rc, result, error
    with torch.no_grad():
      for seed in (17,29):
        oldmodels, oldhash, normalizer = load_models(seed,13312,26624)
        chart = views.load_chart(normalizer); scale = chart['coarse_logrho']
        for case in prep['cases']:
            cap = case['cap']; directory = ROOT/cap
            if sha(directory/'conditions.npz') != case['condition_sha256']:
                raise ValueError('Conditions changed')
            with np.load(directory/'conditions.npz') as f:
                arrays = {k:f[k] for k in f.files}
            with np.load(directory/'positions.npz') as f:
                position = f['position'].T
            tensor = lambda x: torch.as_tensor(x,device='cuda',dtype=torch.float32)[None]
            wide_obs = tensor(views.scale_channels(arrays['wide'],chart['wide']))
            base = FieldCondition(tensor(views.scale_channels(arrays['old_joint'],chart['joint'])),
                                  wide_obs,tensor(np.zeros(3,dtype='f4')),'wide')
            coarse = []; noise = []; coarse_seconds = []
            for draw in range(4):
                tick = time.monotonic()
                z = oldsample(oldmodels['coarse'],base,'cfm',64,[int(c.digest(['loa-page-pilot',cap,seed,draw,'coarse'])[:15],16)])
                torch.cuda.synchronize(); coarse_seconds.append(time.monotonic()-tick)
                coarse.append(z)
                rng = torch.Generator(device='cuda').manual_seed(int(c.digest(['loa-page-pilot',cap,seed,draw,'fine'])[:15],16))
                noise.append(torch.randn((1,1,*SHAPE),device='cuda',generator=rng))
            for architecture,prediction in [('legacy','v'),('unet','x'),('transformer','v')]:
                if time.monotonic()-start > 6300:
                    raise TimeoutError('Pilot reached its conservative work deadline; no automatic retry')
                modelid = f'{architecture}_{prediction}_{seed}'
                if architecture == 'legacy':
                    model = oldmodels['fine']; checkpoint = oldhash['fine']
                else:
                    path = FITROOT/'fits'/modelid/'CHECKPOINT_009984.pt'
                    state = torch.load(path,map_location='cpu',weights_only=False)
                    if state['binding']['sources'] != sources() or state['binding']['normalizer'] != normalizer or state['step'] != 9984:
                        raise ValueError('Expanded checkpoint binding mismatch')
                    model = FineModel(architecture,prediction).cuda().eval(); model.load_state_dict(state['ema'])
                    checkpoint = sha(path); del state
                output = directory/modelid; output.mkdir(exist_ok=False)
                timings = []; checks = []; core_draws = []; tidal_draws = []; galaxy_eigen = []; galaxy_rho = []
                torch.cuda.reset_peak_memory_stats()
                for draw,z in enumerate(coarse):
                    tick = time.monotonic()
                    wide = np.exp(z[0,0].double().cpu().numpy()*scale['std'][0]+scale['mean'][0])
                    if architecture == 'legacy':
                        cut = tuple(slice(*v) for v in op.layout()['joint_coarse_crop_in_wide'])
                        condition = replace(base,region='joint',coarse_joint=replicate(z[(slice(None),slice(None),*cut)]),coarse_wide=z,coarse_source='sampled')
                        # Legacy public sampler owns addressed noise. Its seed is
                        # retained, but noise is not claimed matched across shapes.
                        latent = oldsample(model,condition,'cfm',64,[int(c.digest(['loa-page-pilot',cap,seed,draw,'legacy-fine'])[:15],16)])
                        mass = torch.as_tensor(wide[cut],device='cuda')[None,None]
                    else:
                        condition = Condition(tensor(views.scale_channels(arrays['joint'],chart['joint'])),wide_obs,z,'sampled')
                        latent = sample(model,condition,noise[draw],128)
                        mass = torch.as_tensor(wide[COARSE_CROP],device='cuda')[None,None]
                    rho = decode_density(mass,latent.double()*chart['fine_residual']['std'][0])[0,0].cpu().numpy()
                    mass_error = float(np.max(abs(op.mean_pool(rho)-mass[0,0].cpu().numpy())/mass[0,0].cpu().numpy()))
                    if not np.isfinite(rho).all() or rho.min() <= 0 or mass_error > 2e-6:
                        raise ValueError('Positive density / coarse mass conservation failure')
                    torch.cuda.synchronize(); sampling_seconds = time.monotonic()-tick
                    rc,tc,trace_error = tensor_core(rho,wide,'old' if architecture == 'legacy' else 'new')
                    gt = np.stack([map_coordinates(tc[...,i],position,order=1,mode='nearest',prefilter=False) for i in range(6)],-1)
                    gr = map_coordinates(rc,position,order=1,mode='nearest',prefilter=False)
                    ge = tidal.eigensystem(gt)
                    np.testing.assert_allclose(ge.sum(-1),gr-1,rtol=2e-7,atol=2e-7)
                    core_draws.append(rc.astype('f4')); tidal_draws.append(tc.astype('f4'))
                    galaxy_eigen.append(ge.astype('f4')); galaxy_rho.append(gr.astype('f4'))
                    checks.append({'mass_relative':mass_error,'trace_maxabs':trace_error})
                    timings.append({'sampling_seconds':sampling_seconds,'including_tides_seconds':time.monotonic()-tick})
                    print(json.dumps({'stage':'draw','cap':cap,'model':modelid,'draw':draw,'seconds':timings[-1]}),flush=True)
                np.savez_compressed(output/'draws.npz',rho=np.asarray(core_draws),tensor=np.asarray(tidal_draws),
                    galaxy_eigen=np.asarray(galaxy_eigen),galaxy_rho=np.asarray(galaxy_rho))
                record = {'passed':True,'model':modelid,'cap':cap,'draws':4,'nfe':128,'checks':checks,'timings':timings,
                    'coarse_seconds':coarse_seconds,'peak_gpu_bytes':torch.cuda.max_memory_allocated(),
                    'checkpoint_sha256':checkpoint,'coarse_checkpoint_sha256':oldhash['coarse'],'normalizer':normalizer,
                    'output_sha256':sha(output/'draws.npz'),'closure':'finite free-space local plus wide, unit density beyond wide; not learned full-cap exterior',
                    'truth_opened':False,'job':os.environ['SLURM_JOB_ID']}
                save_json(output/'COMPLETE.json',record); all_records.append(record)
                if architecture != 'legacy':
                    del model; torch.cuda.empty_cache()
            del coarse, noise
        del oldmodels; torch.cuda.empty_cache()
    save_json(ROOT/'INFERENCE_COMPLETE.json',{'passed':True,'records':all_records,'seconds':time.monotonic()-start,
        'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'runner_sha256':sha(__file__),
        'scope':'Two finite real-Loa domains; six frozen fine fits; four draws each; exploratory, not calibrated or full survey'})
    print('GPU_PILOT_COMPLETE',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','infer'])
    action = parser.parse_args().action
    (prepare if action == 'prepare' else infer)()
