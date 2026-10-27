"""Observed-Loa six-fit runner. No scheduler submission, truth access or training.

prepare defaults to the 16-tile geometry-selected preflight. Inference is an
explicit half-open draw range; complete, source-bound draw bundles are reused.
Independent tile draws are NOT a joint survey posterior. Production requires
separate human authorization; this module never submits or retries a job.
"""
import os
for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_name] = '4'
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
from pathlib import Path
import resource
import socket
import sys
import tempfile
import time

ILLUSTRIS = Path('/global/u2/d/dkololgi/TNG/Illustris')
sys.path.insert(0, str(ILLUSTRIS))
ATLAS = Path('/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1')
VAC = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')
PHOT = Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
DEFAULT_ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_preflight_20261016_v2')
SEEDS = (17, 29)
FAMILIES = (('legacy', 'v'), ('unet', 'x'), ('transformer', 'v'))
SHELLS = (.15, .25, .35, .45, .55)
CLOSURE = 'finite free-space local plus wide; unit density beyond wide; independent tiles, not a joint survey posterior'


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 << 20), b''):
            result.update(block)
    return result.hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def require_compute():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Payload work requires an allocated compute node')


def atomic_json(path, value):
    """Publish new immutable metadata; never replace an existing receipt."""
    path = Path(path)
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    with os.fdopen(fd, 'w') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.flush()
        os.fsync(stream.fileno())
    os.link(temporary, path)  # atomic no-clobber publication
    os.unlink(temporary)


@contextmanager
def lock(path, shared=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('a') as stream:
        fcntl.flock(stream, (fcntl.LOCK_SH if shared else fcntl.LOCK_EX) | fcntl.LOCK_NB)
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def verified_bundle(path, binding):
    if not path.exists():
        return None
    receipt = read(path/'COMPLETE.json')
    if receipt.get('binding') != binding or not receipt.get('passed'):
        raise ValueError('Bundle binding mismatch: ' + str(path))
    if sha(path/'arrays.npz') != receipt['sha256']:
        raise ValueError('Bundle checksum mismatch: ' + str(path))
    return receipt


def commit_bundle(path, arrays, metadata, binding):
    """One atomic directory rename commits arrays and receipt together.

    Call under the enclosing tile/index lock. Incomplete hidden temporary
    directories are never treated as output or automatically deleted.
    """
    import numpy as np
    if path.exists():
        raise FileExistsError(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix='.' + path.name + '.', dir=path.parent))
    with (temporary/'arrays.npz').open('xb') as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    receipt = dict(metadata, binding=binding, passed=True, sha256=sha(temporary/'arrays.npz'),
                   bytes=(temporary/'arrays.npz').stat().st_size)
    atomic_json(temporary/'COMPLETE.json', receipt)
    os.rename(temporary, path)
    return receipt


def source_binding():
    """Bind the runner plus local scientific Python implementation dependencies."""
    files = [Path(__file__)]
    # Include the small source tree, not data/products. This also covers indirect
    # imports that the checkpoint's narrower training-source binding omits.
    files += sorted((ILLUSTRIS/'workflows/sbi').glob('e2e_*.py'))
    files += [ILLUSTRIS/'workflows/abacus_tweb'/name for name in
              ('p3a_build_canonical_fields.py', 'p3br_build_random_response.py')]
    return {str(path): sha(path) for path in files}


def checkpoint_binding():
    import torch
    from workflows.sbi.e2e_cfm_pilot_evaluate import TRAIN_ROOT
    from workflows.sbi.e2e_wholefield_pilot_data import ROOT as FITROOT
    from workflows.sbi.e2e_wholefield_pilot_train import sources
    result = {}
    for seed in SEEDS:
        for stage, step in [('coarse', 13312), ('fine', 26624)]:
            path = TRAIN_ROOT/f'{stage}_{seed}'/f'CHECKPOINT_{step:06d}.pt'
            state = torch.load(path, map_location='cpu', weights_only=False)
            if state['step'] != step or state['binding']['seed'] != seed or state['binding']['confirmation_access']:
                raise ValueError('Legacy checkpoint identity mismatch')
            result[f'{stage}_{seed}'] = {'path':str(path), 'sha256':sha(path)}
            del state
        for architecture, prediction in FAMILIES[1:]:
            modelid = f'{architecture}_{prediction}_{seed}'
            path = FITROOT/'fits'/modelid/'CHECKPOINT_009984.pt'
            state = torch.load(path, map_location='cpu', weights_only=False)
            if state['step'] != 9984 or state['binding']['sources'] != sources():
                raise ValueError('Expanded checkpoint source mismatch')
            result[modelid] = {'path':str(path), 'sha256':sha(path)}
            del state
    return result


def shell_index(z):
    if not SHELLS[0] <= z <= SHELLS[-1]:
        raise ValueError('Redshift outside supported range')
    return min(3, int((z - .15 + 1e-12) // .1))


def select_preflight(cases):
    """Sixteen distinct tiles, selected before inspecting fields/properties.

    Shell is the median observed galaxy redshift. Edge/interior are lower/upper
    quartiles of median catalogue BOUNDARY_MPC within each cap/shell. These are
    relative strata, not a guarantee of pristine interior support.
    """
    chosen = []
    for cap in ('NGC', 'SGC'):
        for shell in range(4):
            candidates = sorted((r for r in cases if r['cap'] == cap and
                                 r['shell'] == shell and r['retained_rows'] > 0),
                                key=lambda r:(r['boundary_median'], r['id']))
            if len(candidates) < 2:
                raise ValueError(f'Insufficient geometric coverage: {cap} shell {shell}')
            for label, quantile in [('edge', .25), ('interior', .75)]:
                index = round(quantile*(len(candidates)-1))
                row = dict(candidates[index], stratum=label)
                chosen.append(row)
    if len({r['id'] for r in chosen}) != 16:
        raise ValueError('Preflight tiles are not unique')
    return chosen


def index_catalogue(root):
    """One complete, unambiguous photometry join reused by all tile preparation."""
    import fitsio
    import h5py
    import numpy as np
    from workflows.sbi.e2e_loa_field_atlas import groups
    root.mkdir(parents=True, exist_ok=True)
    marker = root/'INDEX.json'
    with lock(root/'index.lock', shared=marker.exists()):
        if marker.exists():
            result = read(marker)
            if result['sources'] != source_binding() or sha(root/'catalogue.h5') != result['catalogue_sha256']:
                raise ValueError('Index source/catalogue drift')
            if sha(ATLAS/'PLAN.json') != result['atlas_plan_sha256']:
                raise ValueError('Atlas plan changed')
            for path, expected in result['response_sources'].items():
                if sha(path) != expected:
                    raise ValueError('Observed response source changed')
            for cache in result['caches'].values():
                path = Path(cache['path'])
                if path.stat().st_size != cache['bytes'] or path.stat().st_mtime_ns != cache['mtime_ns']:
                    raise ValueError('Observation cache changed after full content verification')
            return result
        start = time.monotonic()
        plan = read(ATLAS/'PLAN.json')
        if sha(ATLAS/'vac_properties.npz') != plan['properties_sha256']:
            raise ValueError('Atlas property source changed')
        if sha(VAC) != plan['property_product']['sha256']:
            raise ValueError('Parent VAC source changed')
        with np.load(ATLAS/'vac_properties.npz') as saved:
            props = {name:saved[name] for name in ('xyz', 'CAP', 'TARGETID', 'Z')}
        columns = ['TARGETID','RA','DEC','Z','SUPPORTED','BOUNDARY_MPC','CIGALE_MATCHED',
                   'CIGALE_AMBIGUOUS','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15',
                   'AGE_CG_15','AV_CG_15','FAGN_CG_15','P_WEB']
        a = fitsio.read(VAC, columns=columns)
        np.testing.assert_array_equal(a['TARGETID'], props['TARGETID'])
        np.testing.assert_array_equal(a['Z'], props['Z'])
        order = np.argsort(a['TARGETID']); ids = a['TARGETID'][order]
        if not np.all(np.diff(ids) > 0):
            raise ValueError('Duplicate parent TARGETID')
        r = np.full(len(a), np.nan); matched = np.zeros(len(a), dtype='i4')
        join_start = time.monotonic()
        with fitsio.FITS(PHOT) as source:
            for first in range(0, source[1].get_nrows(), 200000):
                b = source[1].read(rows=np.arange(first, min(first+200000, source[1].get_nrows())),
                    columns=['TARGETID','RA','DEC','Z_not4clus','FLUX_R','MW_TRANSMISSION_R','PHOTSYS'])
                at = np.searchsorted(ids, b['TARGETID'])
                keep = (at < len(ids)) & (ids[np.minimum(at, len(ids)-1)] == b['TARGETID'])
                ii = order[at[keep]]; b = b[keep]
                for left, right in [('RA','RA'),('DEC','DEC'),('Z_not4clus','Z')]:
                    np.testing.assert_array_equal(a[right][ii], b[left])
                with np.errstate(divide='ignore', invalid='ignore'):
                    flux = b['FLUX_R']/b['MW_TRANSMISSION_R']
                ps = np.char.strip(b['PHOTSYS'].astype('U'))
                valid = np.isfinite(flux) & (flux > 0) & (b['MW_TRANSMISSION_R'] > 0) & np.isin(ps, ['N','S'])
                r[ii[valid]] = 22.5-2.5*np.log10(flux[valid])-np.where(ps[valid] == 'N', .04, 0)
                np.add.at(matched, ii, 1)  # count duplicates even within a chunk
        if not np.all(matched == 1):
            raise ValueError('Photometry join is missing or ambiguous')
        join_seconds = time.monotonic()-join_start
        grouped = groups(props, plan); cases = []
        ordered_rows = []
        for source in sorted(plan['tiles'], key=lambda v:v['id']):
            case = dict(source); rows = grouped[case['id']]
            if len(rows) != case['rows']:
                raise ValueError('Tile row count mismatch')
            grid = plan['grids'][case['cap']]
            origin = np.asarray(grid['origin_mpc_h'])+(np.asarray(case['midpoint'])-[32,16,16])*3.383
            position = (props['xyz'][rows]-origin)/6.766-.5
            inside = ((position >= 0) & (position <= [31,15,15])).all(1)
            retained = rows[inside]
            boundary = a['BOUNDARY_MPC'][rows]
            if not np.isfinite(boundary).all():
                raise ValueError('Nonfinite geometric boundary distance')
            case.update(core_origin_mpc_h=origin.tolist(), retained_rows=int(inside.sum()),
                        excluded_interpolation_edge=int((~inside).sum()),
                        z_median=float(np.median(a['Z'][rows])),
                        z_range=[float(a['Z'][rows].min()),float(a['Z'][rows].max())],
                        boundary_median=float(np.median(boundary)),
                        row_start=sum(len(v) for v in ordered_rows), row_count=len(retained))
            case['shell'] = shell_index(case['z_median'])
            ordered_rows.append(retained); cases.append(case)
        order_rows = np.concatenate(ordered_rows)
        if len(np.unique(order_rows)) != len(order_rows):
            raise ValueError('Duplicate tile ownership')
        temporary = root/'catalogue.h5.partial'
        with h5py.File(temporary, 'x') as output:
            output.create_dataset('xyz', data=props['xyz'][order_rows], compression='lzf')
            output.create_dataset('r_bright', data=r[order_rows], compression='lzf')
            for name in columns:
                value = a[name][order_rows]
                if value.dtype.kind == 'U':
                    value = value.astype('S')
                output.create_dataset(name, data=value, compression='lzf')
        os.link(temporary, root/'catalogue.h5'); temporary.unlink()
        checkpoints = checkpoint_binding()
        caches = {}
        for cap in ('NGC','SGC'):
            receipt = read(ATLAS/f'{cap}_CACHE.json')
            record = receipt['outputs'][0]
            if sha(record['path']) != record['sha256']:
                raise ValueError('Observation cache changed')
            path = Path(record['path'])
            caches[cap] = {'receipt_sha256':sha(ATLAS/f'{cap}_CACHE.json'), 'data_sha256':record['sha256'],
                           'path':str(path),'bytes':path.stat().st_size,'mtime_ns':path.stat().st_mtime_ns}
        # Bind ancillary observed response inputs, not only the cached counts.
        from workflows.sbi.e2e_cfm_field_gallery import LOA
        from workflows.sbi import e2e_coupled_contract as contract
        response_paths = [LOA/'random_angular.npz', LOA/'boundary_angles.npz',
                          Path(contract.config()['observation']['fixed_selection_manifest'])]
        result = dict(schema='loa-survey-runner-v1', passed=True, cases=cases,
            preflight=select_preflight(cases), grids=plan['grids'], parent_counts=plan['counts'],
            sources=source_binding(), checkpoints=checkpoints, caches=caches,
            response_sources={str(p):sha(p) for p in response_paths},
            atlas_plan_sha256=sha(ATLAS/'PLAN.json'), vac_sha256=plan['property_product']['sha256'],
            photometry={'path':str(PHOT),'bytes':PHOT.stat().st_size,'mtime_ns':PHOT.stat().st_mtime_ns},
            catalogue_sha256=sha(root/'catalogue.h5'), join_seconds=join_seconds,
            seconds=time.monotonic()-start, peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            job=os.environ['SLURM_JOB_ID'], host=socket.gethostname(), truth_opened=False,
            closure=CLOSURE, shell_assignment='median observed redshift; shell edges .15,.25,.35,.45,.55',
            edge_interior='25th/75th percentile median catalogue BOUNDARY_MPC within cap/shell')
        atomic_json(marker, result)
        print(json.dumps({'stage':'index_complete','seconds':result['seconds'],'join_seconds':join_seconds,
                          'retained':len(order_rows),'preflight':[r['id'] for r in result['preflight']]}), flush=True)
        return result


def selected_cases(index, args):
    cases = index['cases'] if args.scope == 'all' else index['preflight']
    if args.tiles:
        wanted = set(args.tiles.split(','))
        cases = [r for r in index['cases'] if r['id'] in wanted]
        if {r['id'] for r in cases} != wanted:
            raise ValueError('Unknown tile IDs')
    if args.cap:
        cases = [r for r in cases if r['cap'] == args.cap]
    if args.shard is not None:
        cases = cases[args.shard*24:(args.shard+1)*24]
    if not cases:
        raise ValueError('Empty tile selection')
    return cases


def prepare(args):
    require_compute()
    import h5py
    import numpy as np
    from scipy.ndimage import map_coordinates
    from workflows.sbi import e2e_loa_field_conditions as lc
    from workflows.sbi import e2e_wholefield_pilot_data as wd
    root = args.root; index = index_catalogue(root)
    source = digest(index); responses = {}
    records = []; started = time.monotonic()
    for case in selected_cases(index, args):
        if time.monotonic()-started > args.seconds:
            raise TimeoutError('Preparation deadline reached; no automatic retry')
        cap = case['cap']; destination = root/'tiles'/case['id']/'prepared'
        with lock(root/'locks'/(case['id']+'.lock')):
            old_receipt = verified_bundle(destination, source)
            if old_receipt:
                records.append(old_receipt); continue
            tick = time.monotonic()
            if cap not in responses:
                responses[cap] = lc.load_response(cap)
            with h5py.File(ATLAS/f'{cap}_observation_cache.h5','r') as cache:
                old = lc.build(index['grids'][cap], case['midpoint'], cap, cache['counts_raw'], cache['wide'], responses[cap])
                expanded = wd.build_fine_observations(index['grids'][cap], case['midpoint'], cap, cache['counts_raw'], responses[cap])
            error = float(np.max(abs(expanded[(slice(None), *wd.OLD_CROP)]-old['joint'])))
            if error > 2e-5:
                raise ValueError('Expanded/legacy observation parity failure')
            with h5py.File(root/'catalogue.h5','r') as catalogue:
                sl = slice(case['row_start'], case['row_start']+case['row_count'])
                arrays = {name:catalogue[name][sl] for name in catalogue}
            position = (arrays['xyz']-np.asarray(case['core_origin_mpc_h']))/6.766-.5
            support = map_coordinates(old['support'], (position+[16,16,16]).T, order=1, mode='nearest', prefilter=False)
            arrays.update(joint=expanded, old_joint=old['joint'], wide=old['wide'], position=position, core_support=support)
            mass = arrays['MASS_CG_15']; galaxy = np.char.strip(arrays['CIGALE_SPECTYPE'].astype('U')) == 'GALAXY'
            eligible = (arrays['SUPPORTED'].astype(bool) & (support >= .5) & arrays['CIGALE_MATCHED'].astype(bool) &
                        ~arrays['CIGALE_AMBIGUOUS'].astype(bool) & galaxy & np.isfinite(mass) &
                        (mass >= 1e9) & (mass < 1e12) & (arrays['r_bright'] < 19.5))
            metadata = dict(tile=case['id'], case=case, seconds=time.monotonic()-tick,
                            central_channel_maxabs=error, property_eligible=int(eligible.sum()),
                            eligible_by_shell=np.histogram(arrays['Z'][eligible], SHELLS)[0].tolist(),
                            retained_by_shell=np.histogram(arrays['Z'], SHELLS)[0].tolist(),
                            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                            job=os.environ['SLURM_JOB_ID'], host=socket.gethostname())
            receipt = commit_bundle(destination, arrays, metadata, source)
            records.append(receipt)
            print(json.dumps({'stage':'prepared', **metadata}), flush=True)
    record = dict(passed=True, index_sha256=sha(root/'INDEX.json'), records=records,
                  seconds=time.monotonic()-started, job=os.environ['SLURM_JOB_ID'])
    atomic_json(root/f'PREPARE_{os.environ["SLURM_JOB_ID"]}_{time.time_ns()}.json', record)


def draw_seed(tile, seed, draw, stage):
    # Tile identity is essential: the two-tile pilot used cap-only addressing.
    # The namespace excludes requested draw count, worker, shard and checkpoint.
    return int(digest(['loa-survey-v1',tile,seed,draw,stage])[:15],16)


def infer(args):
    require_compute()
    import numpy as np
    import torch
    from dataclasses import replace
    from scipy.ndimage import map_coordinates
    from workflows.sbi import e2e_coupled_views as views
    from workflows.sbi import e2e_coupled_operators as op
    from workflows.sbi import e2e_footprint_tidal_operator as tidal
    from workflows.sbi.e2e_cfm_pilot_evaluate import load_models
    from workflows.sbi.e2e_coupled_benchmark_models import FieldCondition, sample as oldsample
    from workflows.sbi.e2e_vdm_context_models import replicate, decode_density
    from workflows.sbi.e2e_wholefield_pilot_models import FineModel, Condition, sample
    from workflows.sbi.e2e_wholefield_pilot_train import setup, sources
    from workflows.sbi.e2e_wholefield_pilot_data import SHAPE, COARSE_CROP
    start = time.monotonic(); root = args.root; index = read(root/'INDEX.json')
    if index['sources'] != source_binding():
        raise ValueError('Scientific source changed since preparation')
    for record in index['checkpoints'].values():
        if sha(record['path']) != record['sha256']:
            raise ValueError('Frozen checkpoint changed')
    setup(17); cases = selected_cases(index, args)
    source = digest(index); kernels = {}; records = []; reused = 0; new = 0
    # Float64 free-space operator exactly as in the two-region pilot.
    for key, shape, offset in [('wide',(192,)*3,(80,88,88)), ('new',SHAPE,(32,24,24)), ('old',(64,48,48),(16,16,16))]:
        k = tidal.kernel(shape, (32,16,16), offset=offset)
        kernels[key] = (k.shape[1:], [torch.fft.rfftn(torch.as_tensor(v,device='cuda')) for v in k])
    cs = (slice(32,64),slice(24,40),slice(24,40)); oldcs = (slice(16,48),slice(16,32),slice(16,32))

    def tensor_core(rho, wide, kind):
        lifted = op.lift(wide)
        cut = (slice(48,144),slice(64,128),slice(64,128)) if kind == 'new' else (slice(64,128),slice(72,120),slice(72,120))
        result = np.zeros((32,16,16,6))
        for value, key in [(lifted-1,'wide'),(rho-lifted[cut],kind)]:
            shape, ks = kernels[key]; spectrum = torch.fft.rfftn(torch.as_tensor(value,device='cuda'),s=shape)
            for i,k in enumerate(ks):
                result[...,i] += torch.fft.irfftn(spectrum*k,s=shape)[:32,:16,:16].cpu().numpy()
        rc = rho[cs if kind == 'new' else oldcs]
        error = float(np.max(abs(result[...,[0,3,5]].sum(-1)-(rc-1))))
        if not np.isfinite(result).all() or error > 2e-8*max(1.,float(abs(rc-1).max())):
            raise ValueError('Core free-space trace failed')
        return rc, result, error

    def deadline():
        if time.monotonic()-start > args.seconds:
            raise TimeoutError('Inference deadline reached; completed draws retained; no automatic retry')

    with torch.no_grad():
      for seed in SEEDS:
        deadline()
        oldmodels, oldhash, normalizer = load_models(seed,13312,26624)
        chart = views.load_chart(normalizer); scale = chart['coarse_logrho']; models = {}
        for architecture, prediction in FAMILIES[1:]:
            modelid = f'{architecture}_{prediction}_{seed}'; ck = index['checkpoints'][modelid]
            state = torch.load(ck['path'],map_location='cpu',weights_only=False)
            if state['binding']['sources'] != sources() or state['binding']['normalizer'] != normalizer:
                raise ValueError('Expanded model/normalizer binding mismatch')
            model = FineModel(architecture,prediction).cuda().eval(); model.load_state_dict(state['ema'])
            models[architecture] = model; del state
        for case in cases:
          with lock(root/'locks'/(case['id']+'.lock')):
            directory = root/'tiles'/case['id']; prep = verified_bundle(directory/'prepared', source)
            if prep is None:
                raise ValueError('Tile not prepared: '+case['id'])
            with np.load(directory/'prepared/arrays.npz') as f:
                arrays = {name:f[name] for name in ('joint','old_joint','wide','position')}
            position = arrays['position'].T
            tensor = lambda x:torch.as_tensor(x,device='cuda',dtype=torch.float32)[None]
            wide_obs = tensor(views.scale_channels(arrays['wide'],chart['wide']))
            base = FieldCondition(tensor(views.scale_channels(arrays['old_joint'],chart['joint'])),
                                  wide_obs,tensor(np.zeros(3,dtype='f4')),'wide')
            new_obs = tensor(views.scale_channels(arrays['joint'],chart['joint']))
            binding = digest(dict(index=source, prepared=prep['sha256'], seed=seed, nfe=128,
                                  noise='loa-survey-v1', normalizer=normalizer, closure=CLOSURE))
            for draw in range(args.first, args.stop):
                deadline(); tick = time.monotonic()
                coarse_path = directory/f'coarse_{seed}'/f'{draw:03d}'
                coarse_record = verified_bundle(coarse_path, binding)
                if coarse_record is None:
                    z = oldsample(oldmodels['coarse'],base,'cfm',64,[draw_seed(case['id'],seed,draw,'coarse')])
                    torch.cuda.synchronize(); coarse_seconds = time.monotonic()-tick
                    if not torch.isfinite(z).all():
                        raise ValueError('Nonfinite coarse sample')
                    coarse_record = commit_bundle(coarse_path, {'z':z.cpu().numpy()},
                        dict(seconds=coarse_seconds, seed=seed, draw=draw, tile=case['id'],
                             checkpoint_sha256=oldhash['coarse'], job=os.environ['SLURM_JOB_ID']), binding)
                else:
                    with np.load(coarse_path/'arrays.npz') as f:
                        z = torch.as_tensor(f['z'],device='cuda')
                wide = np.exp(z[0,0].double().cpu().numpy()*scale['std'][0]+scale['mean'][0])
                if not np.isfinite(wide).all() or np.any(wide <= 0):
                    raise ValueError('Invalid coarse density')
                rng = torch.Generator(device='cuda').manual_seed(draw_seed(case['id'],seed,draw,'fine'))
                noise = torch.randn((1,1,*SHAPE),device='cuda',generator=rng)
                for architecture,prediction in FAMILIES:
                    deadline(); modelid = f'{architecture}_{prediction}_{seed}'
                    output = directory/modelid/f'{draw:03d}'
                    draw_binding = digest([binding,modelid,draw,coarse_record['sha256']])
                    previous = verified_bundle(output, draw_binding)
                    if previous:
                        records.append(previous); reused += 1; continue
                    tick = time.monotonic(); torch.cuda.reset_peak_memory_stats()
                    if architecture == 'legacy':
                        cut = tuple(slice(*v) for v in op.layout()['joint_coarse_crop_in_wide'])
                        condition = replace(base,region='joint',coarse_joint=replicate(z[(slice(None),slice(None),*cut)]),coarse_wide=z,coarse_source='sampled')
                        latent = oldsample(oldmodels['fine'],condition,'cfm',64,[draw_seed(case['id'],seed,draw,'legacy-fine')])
                        mass = torch.as_tensor(wide[cut],device='cuda')[None,None]; checkpoint = oldhash['fine']
                    else:
                        condition = Condition(new_obs,wide_obs,z,'sampled')
                        latent = sample(models[architecture],condition,noise,128)
                        mass = torch.as_tensor(wide[COARSE_CROP],device='cuda')[None,None]
                        checkpoint = index['checkpoints'][modelid]['sha256']
                    rho = decode_density(mass,latent.double()*chart['fine_residual']['std'][0])[0,0].cpu().numpy()
                    mass_error = float(np.max(abs(op.mean_pool(rho)-mass[0,0].cpu().numpy())/mass[0,0].cpu().numpy()))
                    if not np.isfinite(rho).all() or rho.min() <= 0 or mass_error > 2e-6:
                        raise ValueError('Density/coarse-mass conservation failure')
                    torch.cuda.synchronize(); sampling_seconds = time.monotonic()-tick
                    rc,tc,trace_error = tensor_core(rho,wide,'old' if architecture == 'legacy' else 'new')
                    gt = np.stack([map_coordinates(tc[...,i],position,order=1,mode='nearest',prefilter=False) for i in range(6)],-1)
                    gr = map_coordinates(rc,position,order=1,mode='nearest',prefilter=False)
                    ge = tidal.eigensystem(gt)
                    np.testing.assert_allclose(ge.sum(-1),gr-1,rtol=2e-7,atol=2e-7)
                    values = dict(rho=rc.astype('f4'),tensor=tc.astype('f4'),galaxy_eigen=ge.astype('f4'),galaxy_rho=gr.astype('f4'))
                    if not all(np.isfinite(v).all() for v in values.values()):
                        raise ValueError('Nonfinite saved float32 values')
                    if np.any(values['rho'] <= 0) or np.any(np.diff(values['galaxy_eigen'],axis=-1) < 0):
                        raise ValueError('Saved density positivity/eigenvalue order failed')
                    metadata = dict(tile=case['id'], model=modelid, seed=seed, draw=draw, nfe=128,
                        sampling_seconds=sampling_seconds, including_tides_seconds=time.monotonic()-tick,
                        mass_relative=mass_error, trace_maxabs=trace_error,
                        peak_gpu_bytes=torch.cuda.max_memory_allocated(), peak_gpu_reserved_bytes=torch.cuda.max_memory_reserved(),
                        peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        coarse_sha256=coarse_record['sha256'], checkpoint_sha256=checkpoint,
                        job=os.environ['SLURM_JOB_ID'], host=socket.gethostname(), truth_opened=False)
                    receipt = commit_bundle(output,values,metadata,draw_binding)
                    receipt['write_seconds'] = time.monotonic()-tick-metadata['including_tides_seconds']
                    records.append(receipt); new += 1
                    print(json.dumps({'stage':'draw_complete', **receipt}),flush=True)
                    if args.max_new and new >= args.max_new:
                        result = dict(passed=True, intentional_stop=True, new=new, reused=reused,
                                      seconds=time.monotonic()-start, records=records)
                        atomic_json(root/f'INFER_{os.environ["SLURM_JOB_ID"]}_{time.time_ns()}.json',result)
                        return
        del oldmodels, models, model; torch.cuda.empty_cache()
    result = dict(passed=True,new=new,reused=reused,records=records,seconds=time.monotonic()-start,
                  first=args.first,stop=args.stop,tiles=[r['id'] for r in cases],job=os.environ['SLURM_JOB_ID'],
                  gpu=torch.cuda.get_device_name(),runner_sha256=sha(__file__),closure=CLOSURE)
    atomic_json(root/f'INFER_{os.environ["SLURM_JOB_ID"]}_{time.time_ns()}.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','infer'])
    parser.add_argument('--root',type=Path,default=DEFAULT_ROOT)
    parser.add_argument('--scope',choices=['preflight','all'],default='preflight')
    parser.add_argument('--tiles',help='Explicit comma-separated atlas tile IDs')
    parser.add_argument('--cap',choices=['NGC','SGC'])
    parser.add_argument('--shard',type=int,help='Zero-based shard of 24 selected tiles')
    parser.add_argument('--first',type=int,default=0)
    parser.add_argument('--stop',type=int,default=4)
    parser.add_argument('--seconds',type=float,default=6300)
    parser.add_argument('--max-new',type=int,default=0,help='Intentional bounded stop for resume testing')
    args = parser.parse_args()
    if not 0 <= args.first < args.stop <= 32 or (args.shard is not None and args.shard < 0) or args.seconds <= 0 or args.max_new < 0:
        parser.error('Invalid draw range/shard/deadline/max-new')
    (prepare if args.action == 'prepare' else infer)(args)


if __name__ == '__main__':
    main()
