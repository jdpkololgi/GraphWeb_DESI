"""Development-only ph006 morphology/truth audit using frozen posterior draws.

Run on an approved compute allocation. No fitting, phase selection or new truth
release. Local slabs supplement (and are never pooled with) the fixed A06 slab.
"""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from compare_offslab import shape, contrast2
from prepare_parent import digest

ROOT = Path('/pscratch/sd/d/dkololgi/abacus/p10_multiphase')
CAND = ROOT / 'p12a_halo48_candidate_20260924_v1'
R = 7 / .6766
CUT = .5847730874520835


def eigen(t):
    a = np.empty_like(t)
    a[..., 0] = t[..., 0]
    a[..., 1] = a[..., 0] + np.logaddexp(0, t[..., 1])
    a[..., 2] = a[..., 1] + np.logaddexp(0, t[..., 2])
    return a


def fraction(y, mask, weight, block):
    mask = mask & np.isfinite(y)
    if not mask.any():
        return {'rows': 0}
    w = weight[mask]; v = y[mask]
    _, b = np.unique(block[mask], axis=0, return_inverse=True)
    sw = np.bincount(b, weights=w); sy = np.bincount(b, weights=w*v)
    rng = np.random.default_rng(171017)
    draws = rng.integers(len(sw), size=(512, len(sw)))
    means = sy[draws].sum(1) / sw[draws].sum(1)
    return {'rows': int(mask.sum()), 'weighted_mean': float(np.average(v, weights=w)),
            'ci95': np.quantile(means, [.025, .975]).tolist(), 'blocks': len(sw),
            'effective_rows': float(w.sum()**2 / (w*w).sum())}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Approved compute allocation required')
    if args.out.exists():
        raise FileExistsError(args.out)
    import torch
    ds = CAND/'dataset'; po = CAND/'posterior'; au = po/'calibration_audit'
    mp = ds/'P12A_DATASET_READY.json'; ready = json.loads(mp.read_text())
    complete = json.loads((po/'P12A_COMPLETE.json').read_text())
    assert ready['pass'] and not ready['sealed_phase_opened'] and not complete['sealed_phase_opened']
    assert complete['dataset_marker_sha256'] == digest(mp)
    sample = ds/'ph006_selection_sample.npz'
    assert digest(sample) == ready['validation']['sha256']
    ck = torch.load(po/'fmpe_estimator.pt', map_location='cpu', weights_only=False)
    assert ck['dataset_marker_sha256'] == digest(mp)
    d = np.load(sample); ix = np.load(au/'evaluation_index.npy')
    draws = np.load(au/'evaluation_samples_scaled.npy', mmap_mode='r')
    assert draws.shape == (len(ix), 512, 3) and len(np.unique(ix)) == len(ix)
    ctx = d['context'][ix]; target = d['parent_node_id'][ix]
    truth = d['truth_eigenvalues'][ix]; w = d['natural_weight'][ix].astype(float)
    np.testing.assert_allclose(eigen(d['theta_softplus'][ix]), truth, atol=2e-5)
    assert np.isin(d['fold'][ix], ready['validation']['evaluation_folds']).all()
    p = np.empty((len(ix), 4))
    for lo in range(0, len(ix), 1000):
        a = eigen(np.asarray(draws[lo:lo+1000], dtype=float)*np.asarray(ck['theta_std']) + np.asarray(ck['theta_mean']))
        c = (a > .2).sum(2)
        p[lo:lo+len(a)] = np.stack([(c == j).mean(1) for j in range(4)], axis=1)
    np.testing.assert_allclose(p.sum(1), 1)
    phase = ROOT/'ph006/p1_canonical'; manifest = json.loads((phase/'manifest.json').read_text())
    marker_path = Path(manifest['source_completion_marker']); marker = json.loads(marker_path.read_text())
    assert digest(marker_path) == manifest['source_completion_marker_sha256']
    assert digest(manifest['points']) == marker['artifacts']['points_sha256']
    assert digest(manifest['index']) == manifest['index_sha256']
    ci = np.load(manifest['index']); points = np.load(manifest['points'], mmap_mode='r')
    np.testing.assert_array_equal(ci['parent_node_id'], np.arange(len(points)))
    contract = json.loads((ROOT/'training_contract/phases/ph006/phase_contract.json').read_text())
    truth_path = Path(contract['target']['path'])
    assert digest(truth_path) == contract['target']['sha256']
    np.testing.assert_allclose(np.load(truth_path, mmap_mode='r')[target], truth, atol=1e-6, rtol=0)
    # Cartesian points are observer-frame Planck18 Mpc; convert to the A06 axes.
    xyz = np.column_stack((-points[:, 1], -points[:, 0], points[:, 2])).astype(float)
    np.testing.assert_array_equal(points[target, 3], d['cap'][ix])
    from astropy.cosmology import Planck18
    np.testing.assert_allclose(np.linalg.norm(xyz[target], axis=1), Planck18.comoving_distance(ctx[:, 3]).value, atol=.15, rtol=0)
    # Match the released DESI tracer radial range; keep invalid-truth context tracers.
    radial = np.linalg.norm(xyz, axis=1); limits = Planck18.comoving_distance([.15, .55]).value
    keep = ci['context'] & (radial >= limits[0]) & (radial < limits[1])
    parent = np.flatnonzero(keep); qxyz = xyz[parent]; tree = cKDTree(qxyz)
    sel = (ctx[:, 3] >= .2) & (ctx[:, 3] < .3)
    ev = ix[sel]; ids = target[sel]; pos = xyz[ids]; p = p[sel]; truth = truth[sel]; w = w[sel]; ctx = ctx[sel]
    block = np.column_stack((d['cap'][ev], d['superblock_id'][ev]))
    t = (truth > .2).sum(1); pred = p.argmax(1); n = len(ids)
    boundary = np.expm1(ctx[:, 6])
    ra = (np.degrees(np.arctan2(-pos[:, 0], -pos[:, 1])) % 360)
    original = (abs(pos[:, 2]) < 10) & (abs(pos[:, 0]) < 180) & (pos[:, 1] >= 800) & (pos[:, 1] <= 1250) & (ra >= 130) & (ra <= 230)
    # A cylinder query is a bounding sphere followed by exact projected/slab cuts.
    data = dict(evaluation_index=ev, parent_node_id=ids, xyz_Mpc=pos, probability=p,
                truth_eigenvalues=truth, natural_weight=w, block=block, original=original)
    report = {'job': os.environ['SLURM_JOB_ID'], 'phase': 'ph006 development; already exposed',
              'targets': n, 'tracers': len(parent), 'original_targets': int(original.sum()), 'rows': []}
    rng = np.random.default_rng(17017)
    for view in ['local', 'original']:
        height = pos[:, 2] if view == 'local' else np.zeros(n)
        allowed = np.ones(n, bool) if view == 'local' else original
        projected = np.full(n, np.nan); projected_count = np.zeros(n, int)
        for i in np.flatnonzero(allowed):
            near = np.asarray(tree.query_ball_point(pos[i], np.hypot(R, 20)), dtype=int)
            delta = qxyz[near] - pos[i]
            use = (parent[near] != ids[i]) & (abs(qxyz[near, 2] - height[i]) < 10) & ((delta[:, :2]**2).sum(1) <= R*R)
            near = near[use]; projected_count[i] = len(near)
            if len(near) >= 8: projected[i] = contrast2(qxyz[near, :2])
        data[view+'_projected'] = projected; data[view+'_projected_count'] = projected_count
        for ir, radius in enumerate([.5*R, R, 1.5*R]):
            full = np.full((n, 3), np.nan); slab = full.copy(); matched = full.copy()
            nf = np.zeros(n, int); ns = nf.copy()
            eligible = allowed & (boundary > max(radius, np.hypot(R,20)))
            for i in np.flatnonzero(eligible):
                near = np.asarray(tree.query_ball_point(pos[i], radius), dtype=int)
                near = near[parent[near] != ids[i]]
                thin = near[abs(qxyz[near, 2] - height[i]) < 10]
                nf[i] = len(near); ns[i] = len(thin)
                full[i], _ = shape(qxyz[near]); slab[i], _ = shape(qxyz[thin])
                if len(thin) >= 6:
                    matched[i] = np.mean([shape(qxyz[rng.choice(near, len(thin), replace=False)])[0] for _ in range(32)], axis=0)
            for name, a in [('full',full),('slab',slab),('matched',matched),('nfull',nf),('nslab',ns),('eligible',eligible)]:
                data[f'{view}_{ir}_{name}'] = a
            line = np.isfinite(projected) & (projected >= CUT)
            cohorts = {'all': allowed, 'true_wall': allowed & (t == 1),
                       'true_filament': allowed & (t == 2),
                       'predicted_wall': allowed & (pred == 1) & (p[:,2] < .2),
                       'line_predicted_wall': allowed & line & (pred == 1) & (p[:,2] < .2),
                       'line_true_wall': allowed & line & (t == 1)}
            for minimum in [6,8,12]:
                for name, cohort in cohorts.items():
                    valid = cohort & eligible & (nf >= minimum) & np.isfinite(full).all(1)
                    pair = valid & (ns >= minimum) & np.isfinite(matched).all(1)
                    rec = {'view':view,'radius_Mpc':radius,'minimum':minimum,'cohort':name,'cohort_rows':int(cohort.sum()),
                           'full_geometry_counts':np.bincount(full[valid].argmax(1),minlength=3).tolist(),
                           'truth_counts':np.bincount(t[valid],minlength=4).tolist(),
                           'projected_line_fraction':fraction(line.astype(float),cohort & eligible & np.isfinite(projected),w,block),
                           'sheet_fraction':fraction((full.argmax(1)==1).astype(float),valid,w,block),
                           'line_fraction':fraction((full.argmax(1)==0).astype(float),valid,w,block),
                           'wall_truth_fraction':fraction((t==1).astype(float),valid,w,block),
                           'filament_truth_fraction':fraction((t==2).astype(float),valid,w,block),
                           'wall_probability':fraction(p[:,1],valid,w,block),
                           'wall_calibration_residual':fraction((t==1)-p[:,1],valid,w,block),
                           'delta_PminusL':fraction((full[:,1]-full[:,0])-(slab[:,1]-slab[:,0]),pair,w,block),
                           'matched_delta_PminusL':fraction((matched[:,1]-matched[:,0])-(slab[:,1]-slab[:,0]),pair,w,block)}
                    report['rows'].append(rec)
            print('finished',view,radius,flush=True)
    args.out.mkdir(parents=True)
    np.savez_compressed(args.out/'geometry.npz', **data)
    sources = [Path(__file__), Path(__file__).with_name('compare_offslab.py'), mp, sample, po/'fmpe_estimator.pt', au/'evaluation_index.npy', au/'evaluation_samples_scaled.npy', phase/'manifest.json', marker_path, truth_path, Path(manifest['points']), Path(manifest['index'])]
    report['sha256'] = {str(s):digest(s) for s in sources}
    report['geometry_sha256'] = digest(args.out/'geometry.npz')
    report['limitations'] = 'One exposed phase, no blind validation or DESI truth. Weighted by saved shell sampling weights. Local slabs are target-centred, fixed global orientation; not exact A06 selection. Original target cohort uses A06 cuts but neighbours use padded parent. 512 cap/superblock bootstrap draws; conditional on frozen model and this phase. Tracer covariance is not the matter tidal operator. No retraining.'
    (args.out/'RESULTS.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')

if __name__ == '__main__':
    main()
