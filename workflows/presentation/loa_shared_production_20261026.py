"""Frozen, observation-only four-draw shared-state cap production.

Uses the exact GPU-tested samplers. Outputs are exploratory density, not a VAC
or calibrated posterior. No overwrite, automatic retry, or automatic resume.
"""
import argparse
import gc
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import tempfile
import time

import numpy as np
import torch

import loa_shared_sixfit_20261026 as s
import loa_shared_field_io_20261026 as io


def check_freeze(root, expected_hash):
    path = root / 'FREEZE.json'
    if io.sha(path) != expected_hash:
        raise ValueError('Frozen launch manifest changed')
    manifest = s.read(path)
    for relative, digest in manifest['snapshot_files'].items():
        if io.sha(root / relative) != digest:
            raise ValueError('Frozen source changed: ' + relative)
    for name, record in manifest['evidence'].items():
        if io.sha(record['path']) != record['sha256']:
            raise ValueError('Engineering evidence changed: ' + name)
    # Scientific imports must be isolated from the mutable working repositories.
    for name, module in list(sys.modules.items()):
        if name == 'workflows' or name.startswith('workflows.') or name.startswith('loa_shared_'):
            filename = getattr(module, '__file__', None)
            if filename and not Path(filename).resolve().is_relative_to(root):
                raise ValueError('Scientific import escaped source snapshot: ' + name)
    return manifest


def checkpoint(path, state, stage, step, binding_hash):
    fd, temporary = tempfile.mkstemp(prefix='.' + path.name + '.', dir=path.parent)
    with os.fdopen(fd, 'wb') as stream:
        torch.save(dict(state=state.detach().cpu(), stage=stage, heun_step=step,
                        nfe=128, binding_sha256=binding_hash,
                        automatic_resume=False), stream)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def support_from_observations(observations, geometry, expanded, coverage):
    if s.views.products.LOCAL_CHANNELS[1] != 'support_random':
        raise ValueError('Support channel changed')
    support = np.zeros(geometry.shape, dtype='f4')
    seen = np.zeros(geometry.shape, dtype=bool)
    for i, item in enumerate(observations):
        value = item[1 if expanded else 0][0, 1].cpu().numpy()
        io.merge_support(support, seen, value, geometry.crop(i)[2:])
    if not np.array_equal(seen, coverage > 0):
        raise ValueError('Support assembly/coverage mismatch')
    return support


def expanded_state(index, family, seed, normalizer):
    prediction = 'x' if family == 'unet' else 'v'
    model_id = f'{family}_{prediction}_{seed}'
    ck = index['checkpoints'][model_id]
    if io.sha(ck['path']) != ck['sha256']:
        raise ValueError('Expanded checkpoint changed')
    state = torch.load(ck['path'], map_location='cpu', weights_only=False)
    binding = state['binding']
    if (state['step'] != 9984 or binding['sources'] != s.training_sources()
            or binding['normalizer'] != normalizer or binding['architecture'] != family
            or binding['prediction'] != prediction or binding['seed'] != seed
            or binding['sealed_access']):
        raise ValueError('Expanded model binding mismatch')
    return state


@torch.no_grad()
def run(args):
    if (not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid')
            or torch.cuda.device_count() != 1):
        raise RuntimeError('One explicitly approved compute-node GPU required')
    root = args.snapshot.resolve()
    manifest = check_freeze(root, args.manifest_sha256)
    config = s.read(root / 'GraphWeb_DESI/configs/loa_shared_production_20261026.json')
    task = f'{args.cap}_seed{args.seed}'
    if (not config['approved'] or task not in config['tasks']
            or config['draws'] != 4 or config['nfe'] != 128):
        raise ValueError('Task outside approved production contract')
    if str(root) != manifest['snapshot_root']:
        raise ValueError('Snapshot relocation requires a new binding')
    out = Path(config['output_root']) / task
    # Existing partial/successful directories are never resumed or overwritten.
    out.mkdir(parents=True, exist_ok=False)
    (out / 'checkpoints').mkdir()
    started = time.monotonic()
    s.settings.setup(False)
    torch.backends.cuda.enable_flash_sdp(True)
    torch.backends.cuda.enable_mem_efficient_sdp(True)
    if torch.cuda.get_device_properties(0).total_memory < 70 * 2**30:
        raise RuntimeError('The measured production requires an 80GB GPU')
    index_path = s.INPUT / 'INDEX.json'
    index = s.read(index_path)
    if io.sha(index_path) != manifest['input_index_sha256']:
        raise ValueError('Prepared input index changed')
    expected_input = s.views.c.digest(index)
    atlas = Path(config['atlas'])
    plan = s.read(atlas / 'PLAN.json')
    if io.sha(atlas / 'PLAN.json') != index['atlas_plan_sha256']:
        raise ValueError('Atlas plan mismatch')
    if index['grids'] != plan['grids'] or index['parent_counts'] != plan['counts']:
        raise ValueError('Prepared observations no longer match parent geometry')
    if io.sha(atlas / 'vac_properties.npz') != plan['properties_sha256']:
        raise ValueError('Parent catalogue changed')
    rows = s.select(index, args.cap, True)
    atlas_rows = sorted((r for r in plan['tiles'] if r['cap'] == args.cap), key=lambda r: r['id'])
    if ([(r['id'], r['midpoint']) for r in rows] != [(r['id'], r['midpoint']) for r in atlas_rows]
            or len(rows) != plan['counts'][args.cap]['tiles']):
        raise ValueError('Incomplete or changed cap tiling')
    mid = [r['midpoint'] for r in rows]
    cg = s.Geometry.make(mid, 'coarse')
    geometries = {expanded: s.geometry(mid, expanded) for expanded in (False, True)}
    masks = {expanded: io.coverage_for(g) for expanded, g in geometries.items()}
    with np.load(atlas / 'vac_properties.npz') as parent:
        points = parent['xyz'][parent['CAP'] == (1 if args.cap == 'NGC' else 0)]
    preflight = {}
    for expanded, fg in geometries.items():
        offset = np.asarray(fg.origin) // 4 - cg.origin
        required = [12, 16, 16] if expanded else [16, 18, 18]
        if not np.array_equal(offset, required):
            raise ValueError('Native coarse/fine crop mismatch')
        origin = np.asarray(plan['grids'][args.cap]['origin_mpc_h']) + np.asarray(fg.origin) * io.CELL
        preflight['expanded' if expanded else 'legacy'] = io.check_coverage(
            masks[expanded], points, origin, plan['counts'][args.cap]['vac_rows'])
    del points
    io.atomic_json(out / 'COVERAGE_PREFLIGHT.json', preflight)
    models, hashes, normalizer = s.evaluate.load_models(args.seed, 13312, 26624)
    for stage, digest in hashes.items():
        if digest != index['checkpoints'][f'{stage}_{args.seed}']['sha256']:
            raise ValueError('Legacy checkpoint changed')
    chart = s.views.load_chart(normalizer)
    for family in ('unet', 'transformer'):
        checked = expanded_state(index, family, args.seed, normalizer)
        del checked
    binding = dict(schema='loa-shared-sixfit-production-v1', cap=args.cap, seed=args.seed,
        draws=4, nfe=128, windows=len(rows), manifest_sha256=args.manifest_sha256,
        index_sha256=manifest['input_index_sha256'], atlas_plan_sha256=index['atlas_plan_sha256'],
        properties_sha256=plan['properties_sha256'], normalizer=normalizer,
        checkpoints=index['checkpoints'], job=os.environ['SLURM_JOB_ID'], host=socket.gethostname(),
        tf32=False, expanded_autocast='bfloat16', fused_attention_allowed=True,
        shared_state=True, posterior_qualified=False, truth_opened=False, retraining=False,
        automatic_resume=False, automatic_retries=False, git=manifest['git'],
        gpu_name=torch.cuda.get_device_name(0), gpu_total_gib=torch.cuda.get_device_properties(0).total_memory / 2**30)
    io.atomic_json(out / 'BINDING.json', binding)
    binding_hash = io.sha(out / 'BINDING.json')

    def progress(label):
        def update(step, state):
            elapsed = time.monotonic() - started
            timed_out = elapsed > config['tasks'][task]['walltime_seconds'] - 600
            if step % 8 == 0 or timed_out:
                checkpoint(out / 'checkpoints' / (label + '.pt'), state, label, step, binding_hash)
            value = dict(stage=label, heun_step=step, heun_steps=64, elapsed_seconds=elapsed)
            io.atomic_json(out / 'PROGRESS.json', value, replace=True)
            print(json.dumps(value), flush=True)
            if timed_out:
                raise TimeoutError('Approved task time exhausted; no automatic retry')
        return update

    tick = time.monotonic()
    obs, inputs = s.load_observations(rows, chart, expected_input)
    obs = [tuple(value.expand(4, *value.shape[1:]) for value in item) for item in obs]
    io.atomic_json(out / 'INPUTS.json', inputs)
    supports = {expanded: support_from_observations(obs, g, expanded, masks[expanded])
                for expanded, g in geometries.items()}
    base = s.coarse_base(obs)
    preparation = time.monotonic() - tick
    torch.cuda.reset_peak_memory_stats()
    tick = time.monotonic()
    z, cw = s.integrate(models['coarse'], 'coarse', base,
        s.addressed_noise(cg, args.cap, args.seed, 'coarse', draws=4), cg, 128, progress('coarse'))
    torch.cuda.synchronize()
    coarse_seconds = time.monotonic() - tick
    coarse_hash = hashlib.sha256(z.cpu().numpy().tobytes()).hexdigest()
    results = []
    for family in s.FAMILIES:
        expanded = family != 'legacy'
        fg = geometries[expanded]
        fine = s.fine_conditions(obs, base, z, cg, expanded)
        model_id = f'legacy_v_{args.seed}'
        parity = None
        if expanded:
            model_id = f'{family}_{"x" if family == "unet" else "v"}_{args.seed}'
            state = expanded_state(index, family, args.seed, normalizer)
            net = s.FineModel(family, 'x' if family == 'unet' else 'v').cuda().eval()
            net.load_state_dict(state['ema'])
            del state
            single = s.geometry([mid[0]], True)
            noise = s.addressed_noise(single, args.cap, args.seed, 'expanded-parity', draws=4)
            direct = s.sample(net, fine[0], noise, 4)
            coupled, _ = s.integrate_expanded(net, [fine[0]], noise, single, 4)
            parity = float((direct - coupled).abs().max())
            if parity >= 1e-6:
                raise ValueError('Native single-window parity failed')
            del direct, coupled, noise
        tick = time.monotonic()
        initial = s.addressed_noise(fg, args.cap, args.seed,
            'expanded-fine' if expanded else 'legacy-fine', draws=4)
        if expanded:
            u, fw = s.integrate_expanded(net, fine, initial, fg, 128, progress(model_id))
        else:
            u, fw = s.integrate(models['fine'], 'fine', fine, initial, fg, 128, progress(model_id))
        torch.cuda.synchronize()
        fine_seconds = time.monotonic() - tick
        if not np.array_equal(fw[0, 0].cpu().numpy(), masks[expanded]):
            raise ValueError('Sampler coverage differs from independently assembled geometry')
        field = io.write_field(out / (model_id + '_FIELD.h5'), z, u, masks[expanded], supports[expanded],
            cg, fg, chart, plan['grids'][args.cap]['origin_mpc_h'], args.cap, model_id)
        result = dict(model=model_id, field=field, fine_seconds=fine_seconds,
            coarse_latent_sha256=coarse_hash, native_parity_maxabs=parity,
            coverage=preflight['expanded' if expanded else 'legacy'], posterior_qualified=False)
        io.atomic_json(out / (model_id + '_COMPLETE.json'), dict(passed=True, **result))
        results.append(result)
        print(json.dumps(result), flush=True)
        del initial, u, fw, fine
        if expanded:
            del net
        gc.collect()
        torch.cuda.empty_cache()
    check_freeze(root, args.manifest_sha256)
    if io.sha(index_path) != manifest['input_index_sha256']:
        raise ValueError('Input index changed during production')
    result = dict(passed=True, cap=args.cap, seed=args.seed, draws=4, nfe=128,
        preparation_seconds=preparation, coarse_seconds=coarse_seconds,
        seconds=time.monotonic() - started, peak_gpu_gib=torch.cuda.max_memory_allocated() / 2**30,
        results=results, binding_sha256=binding_hash, posterior_qualified=False,
        verification_scope='Geometry, positive density, native sampler parity and mass only; not calibration')
    io.atomic_json(out / 'COMPLETE.json', result)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', type=Path, required=True)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--cap', choices=('NGC', 'SGC'), required=True)
    parser.add_argument('--seed', choices=(17, 29), type=int, required=True)
    run(parser.parse_args())
