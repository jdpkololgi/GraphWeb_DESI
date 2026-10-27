"""Read-only scientific diagnostics of the frozen four-draw Loa preview.

Writes only a new review JSON, never campaign gates or inference products.
Class fractions count galaxy-draws, not independent galaxies or truths.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
import os
from pathlib import Path
import socket
import tarfile
import time

import numpy as np

ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1')
MODELS = [f'{f}_{s}' for s in (17, 29) for f in ('legacy_v', 'unet_x', 'transformer_v')]
PAIRS = [(0, 1), (0, 2), (3, 4), (3, 5), (0, 3), (1, 4), (2, 5)]
SHELLS = np.array([.15, .25, .35, .45, .55])
# Production comparisons occur in float32. Preserve that exact threshold after
# promoting arrays to float64 for reductions (especially eigenvalues on 0.2).
THRESHOLD = float(np.float32(.2))


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def masks(prep):
    z = prep['Z']
    shell = np.minimum(np.searchsorted(SHELLS, z, side='right') - 1, 3)
    if not np.all((z >= .15) & (z <= .55)):
        raise ValueError('Unexpected redshift')
    support = prep['core_support']
    if not np.isfinite(support).all():
        raise ValueError('Nonfinite support')
    supported = prep['SUPPORTED'].astype(bool)
    eligible = (supported & (support >= .5) & prep['CIGALE_MATCHED'].astype(bool)
                & ~prep['CIGALE_AMBIGUOUS'].astype(bool)
                & (np.char.strip(prep['CIGALE_SPECTYPE'].astype('U')) == 'GALAXY')
                & np.isfinite(prep['MASS_CG_15']) & (prep['MASS_CG_15'] >= 1e9)
                & (prep['MASS_CG_15'] < 1e12) & (prep['r_bright'] < 19.5))
    bands = {'unsupported': ~supported, 'support_lt_0.5': supported & (support < .5),
             'support_0.5_0.8': supported & (support >= .5) & (support < .8),
             'support_ge_0.8': supported & (support >= .8)}
    result = {'all': np.ones(len(z), dtype=bool), 'eligible': eligible}
    result.update(bands)
    for i in range(4):
        result[f'z{i}'] = shell == i
        result[f'z{i}_eligible'] = (shell == i) & eligible
        for name, band in bands.items():
            result[f'z{i}_{name}'] = (shell == i) & band
    return result


def summarize(eigen, rho, mask):
    n = int(mask.sum())
    e, r = eigen[:, mask], rho[:, mask]
    labels = (e > THRESHOLD).sum(axis=-1)
    # Spatial contrast is not posterior uncertainty; record both separately.
    spread = (e[..., 2] - e[..., 0]) / np.sqrt(2)
    probs = np.stack([(labels == k).mean(axis=0) for k in range(4)], axis=-1)
    return {
        'n': n, 'class_counts': np.bincount(labels.ravel(), minlength=4).tolist(),
        'eigen_sum': e.sum(axis=(0, 1)).tolist(), 'rho_sum': float(r.sum()),
        'tidal_spread_sum': float(spread.sum()),
        'posterior_eigen_variance_sum': np.var(e, axis=0, ddof=1).sum(axis=0).tolist(),
        'class_ambiguity_sum': float((1 - probs.max(axis=-1)).sum()),
    }


def pair_summary(e1, e2, mask):
    a, b = e1[:, mask], e2[:, mask]
    delta = a.mean(axis=0) - b.mean(axis=0)
    la, lb = (a > THRESHOLD).sum(axis=-1), (b > THRESHOLD).sum(axis=-1)
    pa = np.stack([(la == k).mean(axis=0) for k in range(4)], axis=-1)
    pb = np.stack([(lb == k).mean(axis=0) for k in range(4)], axis=-1)
    return {'n': int(mask.sum()), 'mean_delta_sum': delta.sum(axis=0).tolist(),
            'mean_delta_sq_sum': (delta**2).sum(axis=0).tolist(),
            'class_probability_tv_sum': float(.5 * np.abs(pa - pb).sum()),
            'matched_draw_class_disagreement': int((la != lb).sum())}


def tile(row, audit):
    name = row['id']
    folder = ROOT / 'packed' / f'{name}_04draws'
    meta = json.loads((folder / 'COMPLETE.json').read_text())
    if sha(folder / 'tile.tar') != audit['archive_sha256']:
        raise ValueError('Archive changed since audit: ' + name)
    with tarfile.open(folder / 'tile.tar', 'r') as tar:
        def load(member, columns):
            content = tar.extractfile(member).read()
            if hashlib.sha256(content).hexdigest() != meta['members'][member]['sha256']:
                raise ValueError('Member checksum: ' + member)
            with np.load(io.BytesIO(content)) as a:
                return {k: a[k] for k in columns}
        prep = load('prepared/arrays.npz', [
            'TARGETID', 'Z', 'core_support', 'SUPPORTED', 'CIGALE_MATCHED',
            'CIGALE_AMBIGUOUS', 'CIGALE_SPECTYPE', 'MASS_CG_15', 'r_bright'])
        groups = masks(prep)
        if len(prep['Z']) != audit['retained'] or groups['eligible'].sum() != audit['eligible']:
            raise ValueError('Population mismatch')
        eigen, density, fields = [], [], []
        for m in MODELS:
            draws = [load(f'{m}/{d:03d}/arrays.npz', ['galaxy_eigen', 'galaxy_rho', 'rho']) for d in range(4)]
            eigen.append(np.stack([x['galaxy_eigen'] for x in draws]).astype(float))
            density.append(np.stack([x['galaxy_rho'] for x in draws]).astype(float))
            fields.append(np.stack([x['rho'] for x in draws]).astype(float))
            original_labels = np.stack([(x['galaxy_eigen'] > .2).sum(axis=-1) for x in draws])
            labels = (eigen[-1] > THRESHOLD).sum(axis=-1)
            np.testing.assert_array_equal(labels, original_labels)
            if np.bincount(labels.ravel(), minlength=4).tolist() != audit['models'][m]['retained_class_draw_counts']:
                raise ValueError('Class reconciliation: ' + name + ' ' + m)
        group_results = {}
        for g, mask in groups.items():
            group_results[g] = {
                'models': {m: summarize(eigen[i], density[i], mask) for i, m in enumerate(MODELS)},
                'pairs': {MODELS[i] + '__' + MODELS[j]: pair_summary(eigen[i], eigen[j], mask) for i, j in PAIRS}}
        # Independently recompute maxima and paired core field differences.
        field = {}
        for i, m in enumerate(MODELS):
            f = fields[i]
            old = np.asarray(audit['models'][m]['draw_core_stats'])
            fresh = np.stack([f.mean(axis=(1, 2, 3)), f.std(axis=(1, 2, 3)),
                              f.min(axis=(1, 2, 3)), f.max(axis=(1, 2, 3))], axis=-1)
            np.testing.assert_allclose(fresh, old, atol=1e-10, rtol=1e-10)
            field[m] = {'max_rho': float(f.max()), 'min_rho': float(f.min()),
                        'spatial_sd_mean': float(f.std(axis=(1, 2, 3)).mean()),
                        'posterior_sd_mean': float(f.std(axis=0, ddof=1).mean())}
        return {'tile': name, 'cap': row['cap'], 'tile_shell': row['shell'],
                'boundary_median': row['boundary_median'],
                'n': len(prep['Z']), 'eligible': int(groups['eligible'].sum()),
                'support_min': float(prep['core_support'].min()) if len(prep['Z']) else None,
                'support_max': float(prep['core_support'].max()) if len(prep['Z']) else None,
                'field': field, 'groups': group_results}


def merge(records):
    result = {}
    for rec in records:
        for group, values in rec['groups'].items():
            # Cap groups remain separate from globally pooled groups.
            for key in (group, rec['cap'] + '_' + group):
                dest = result.setdefault(key, {'models': {}, 'pairs': {}})
                for kind in ('models', 'pairs'):
                    for name, vals in values[kind].items():
                        target = dest[kind].setdefault(name, {})
                        for metric, val in vals.items():
                            if isinstance(val, list):
                                target[metric] = [x + y for x, y in zip(target.get(metric, [0] * len(val)), val)]
                            else:
                                target[metric] = target.get(metric, 0) + val
    return result


def selftest():
    e = np.array([[[0., .1, .3]], [[0., .1, .3]], [[0., .3, .4]], [[0., .3, .4]]])
    x = summarize(e, 1 + e.sum(axis=-1), np.array([True]))
    assert x['class_counts'] == [0, 2, 2, 0]
    assert x['class_ambiguity_sum'] == .5
    boundary = np.array([.2], dtype=np.float32)
    assert not (boundary.astype(float) > THRESHOLD).any()
    assert pair_summary(e, e, np.array([True]))['class_probability_tv_sum'] == 0
    assert summarize(e, e.sum(axis=-1), np.array([False]))['n'] == 0
    p = dict(Z=np.array([.15, .25, .35, .45, .55]), core_support=np.array([.49, .5, .8, 1., 1.]),
             SUPPORTED=np.array([1, 1, 1, 1, 0]), CIGALE_MATCHED=np.ones(5),
             CIGALE_AMBIGUOUS=np.zeros(5), CIGALE_SPECTYPE=np.array(['GALAXY'] * 5),
             MASS_CG_15=np.ones(5) * 1e10, r_bright=np.ones(5) * 18)
    g = masks(p)
    assert [g[f'z{i}'].sum() for i in range(4)] == [1, 1, 1, 2]
    assert g['eligible'].sum() == 3
    assert sum(g[k].sum() for k in ('unsupported', 'support_lt_0.5', 'support_0.5_0.8', 'support_ge_0.8')) == 5
    print('selftest passed', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    if not socket.gethostname().startswith('nid') or not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Archive analysis requires an allocated compute node')
    selftest()
    start = time.monotonic()
    audit = json.loads((ROOT / 'PREVIEW_AUDIT.json').read_text())
    cases = json.loads((ROOT / 'INDEX.json').read_text())['cases']
    by_tile = {x['tile']: x for x in audit['per_tile']}
    selected = cases[:args.limit] if args.limit else cases
    results = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        for n, result in enumerate(pool.map(lambda c: tile(c, by_tile[c['id']]), selected), 1):
            results.append(result)
            if n % 100 == 0:
                print(f'{n}/{len(selected)} tiles {time.monotonic()-start:.1f}s', flush=True)
    combined = merge(results)
    if not args.limit:
        for model in MODELS:
            v = combined['all']['models'][model]
            assert v['n'] == audit['retained']
            assert v['class_counts'] == audit['models'][model]['retained_class_draw_counts']
            assert combined['eligible']['models'][model]['n'] == audit['property_eligible']
    # Keep bounded per-tile diagnostics, not galaxy identities.
    per_tile = [{**{k: v for k, v in x.items() if k != 'groups'},
                 'all': x['groups']['all'], 'eligible_summary': x['groups']['eligible']} for x in results]
    output = {'passed': True, 'extension_authorized': False, 'tiles': len(results),
              'summary': combined, 'per_tile': per_tile,
              'audit_sha256': sha(ROOT / 'PREVIEW_AUDIT.json'), 'source_sha256': sha(__file__),
              'job': os.environ['SLURM_JOB_ID'], 'seconds': time.monotonic() - start,
              'notes': ['Observed-data diagnostics, not calibration.',
                        'Galaxy-draw fractions; no independence assumption or inferential error bars.',
                        'Four-draw posterior variances and probability differences are noisy.',
                        'Seed comparisons include different coarse and fine noise.',
                        'Support strata are diagnostic, not new adopted selection cuts.']}
    with args.output.open('x') as f:
        json.dump(output, f, allow_nan=False)
    print(json.dumps({k: output[k] for k in ('passed', 'tiles', 'seconds', 'job')}), flush=True)


if __name__ == '__main__':
    main()
