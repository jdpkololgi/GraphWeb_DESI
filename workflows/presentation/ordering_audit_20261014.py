"""Bounded identity audit of saved products; no inference or full SED scan."""
import hashlib
import json
from pathlib import Path

import fitsio
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18

ROOT = Path('/global/homes/d/dkololgi/GraphWeb_DESI')
OUT = Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010/legacy_gat_comparison_20261014')
BASE = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs')
RUN = BASE / 'p12a_loa_full_20260925_v1'
INPUT = BASE / 'p12a_loa_canary_20260925_v1'


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def equal(a, b):
    return np.array_equal(a, b, equal_nan=True)


def main():
    report = {'scope': 'Full ID census; sampled saved prediction/coordinate/property chain. No model replay.'}
    manifest = read(RUN / 'SOURCE_MANIFEST.json')
    checked = []
    for name, expected in manifest['files'].items():
        assert sha(RUN / 'source' / name) == expected, name
        if name.endswith(('production.py', 'production_inference.py', 'loa_trial.py', 'p12a_observation_patch.py')):
            assert sha(ROOT / name) == expected, name
            checked.append(name)
    report['frozen_sources_verified'] = checked
    science = Path(manifest['science_source'])
    frozen = read(science.parent / 'RUN_MANIFEST.json')['source_sha256']
    science_files = ['workflows/abacus_tweb/p8_train_unet_patch.py',
                     'workflows/sbi/p12_train_base_response_fmpe.py']
    for name in science_files:
        assert sha(science / name) == frozen[name]
        assert sha(ROOT.parent / 'TNG/Illustris' / name) == frozen[name]
    report['frozen_science_sources_verified'] = science_files

    # Reconstruct the historical catalogue reorder independently of its loader.
    old = pd.read_pickle(ROOT / 'DESI_BGS_PRERELEASE_VAC.pkl')
    d = fitsio.read('/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/loa-combined-lowz.fits',
                     columns=['TARGETID', 'RA', 'DEC', 'Z', 'LOGMSTAR', 'ZWARN', 'DELTACHI2', 'BGS_TARGET'])
    d = d[(d['ZWARN'] == 0) & (d['DELTACHI2'] >= 25) & (d['LOGMSTAR'] >= 9) & (d['BGS_TARGET'] != 0)]
    north = SkyCoord(d['RA'], d['DEC'], unit='deg').galactic.b.deg > 0
    ordered = np.concatenate([d[north], d[~north]])
    fields = ['TARGETID', 'RA', 'DEC', 'Z', 'LOGMSTAR']
    checks = {n: equal(old[n].to_numpy(), ordered[n]) for n in fields}
    assert all(checks.values())
    report['legacy'] = dict(rows=len(old), north_rows=int(north.sum()), south_rows=int((~north).sum()),
                            north_south_order_exact=checks,
                            original_order_id_matches=int(np.count_nonzero(old.TARGETID.to_numpy() == d['TARGETID'])))
    del d, ordered, old

    ready = read(INPUT / 'INPUTS_READY.json')
    inp = np.load(INPUT / 'catalogue.npz')
    ids = inp['targetid']; shell = inp['shell']; order = np.argsort(ids); sid = ids[order]
    complete = read(RUN / 'FULL_VAC_COMPLETE.json')
    vac = Path(complete['catalogue'])
    v_ids = fitsio.read(vac, columns=['TARGETID'])['TARGETID']
    assert np.array_equal(v_ids, np.sort(ids[shell >= 0]))
    assert np.all(np.diff(v_ids) > 0)
    report['full_vac_identity_census'] = dict(rows=len(v_ids), unique_exact_input_ids=True)

    # Choose one small saved core in each cap/redshift shell, using core centres.
    selection = read(Path(ready['selection_path']))
    candidates = {}
    for item in read(RUN / 'DRAW_INDEX.json')['shards']:
        if not 32 <= item['rows'] <= 256:
            continue
        p = Path(item['path']); marker = p.with_name(p.name.replace('_draws.npz', '_COMPLETE.json'))
        case = read(marker)['case']; cap = case['cap']; grid = ready['grids']['NGC' if cap else 'SGC']
        centre = np.asarray(grid['origin_mpc']) + (np.asarray(case['core_start']) + np.asarray(case['core_stop'])) * 2.5
        z = np.interp(np.linalg.norm(centre), selection['cosmology']['radius_grid_mpc'], selection['cosmology']['redshift_grid'])
        sh = int(np.searchsorted([.15, .25, .35, .45, .55], z, side='right') - 1)
        if 0 <= sh <= 3 and (cap, sh) not in candidates:
            candidates[cap, sh] = item
        if len(candidates) == 8:
            break
    assert len(candidates) == 8
    xyz = inp['xyz']; cap_array = inp['cap']; ra = inp['ra']; dec = inp['dec']; redshift = inp['z']
    records = []; sampled = []
    for (cap, sh), item in sorted(candidates.items()):
        path = Path(item['path']); assert sha(path) == item['sha256']
        marker = read(path.with_name(path.name.replace('_draws.npz', '_COMPLETE.json')))
        case = marker['case']; grid = ready['grids']['NGC' if cap else 'SGC']
        table = fitsio.read(path.with_name(path.name.replace('_draws.npz', '.fits')))
        with np.load(path) as f:
            di = f['TARGETID']; draws = f['eigenvalue_draws']
        assert np.array_equal(di, table['TARGETID'])
        ix = order[np.searchsorted(sid, di)]; assert np.array_equal(ids[ix], di)
        assert np.all(np.diff(ix) > 0)  # flatnonzero source-row order
        for name, values in [('RA', ra), ('DEC', dec), ('Z', redshift), ('CAP', cap_array)]:
            assert equal(table[name], values[ix]), name
        voxel = np.floor((xyz[ix] - np.asarray(grid['origin_mpc'])) / 5).astype('i8')
        assert np.all((voxel >= case['core_start']) & (voxel < case['core_stop']))
        assert np.all(cap_array[ix] == cap) and np.all(shell[ix] >= 0)
        yes = table['SUPPORTED']; dd = draws[yes]
        assert np.isnan(draws[~yes]).all()
        probs = np.stack([((dd > .2).sum(-1) == j).mean(1) for j in range(4)], axis=1).astype('f4')
        assert equal(probs, table['P_WEB'][yes])
        assert equal(dd.mean(1), table['EIGENVALUE_MEAN'][yes])
        vi = np.searchsorted(v_ids, di); merged = fitsio.read(vac, rows=vi)
        for name in table.dtype.names:
            if name != 'QUALITY':  # merged VAC adds the documented z>=.35 bit
                assert equal(table[name], merged[name]), name
        assert equal(merged['QUALITY'], table['QUALITY'] | np.where(table['Z'] >= .35, 64, 0).astype('u2'))
        sampled.extend(ix.tolist())
        records.append(dict(cap=cap, centre_shell=sh, core_id=item['core_id'], rows=len(di),
                            actual_shells=np.unique(shell[ix]).tolist(), supported=int(yes.sum())))
    # Independent sky-to-Cartesian reconstruction checks what positions IDs carry.
    ix = np.asarray(sampled); zg = np.linspace(0, .85, 17001)
    radius = np.interp(redshift[ix], zg, Planck18.comoving_distance(zg).value)
    rr = np.deg2rad(ra[ix]); dd = np.deg2rad(dec[ix])
    rebuilt = radius[:, None] * np.column_stack([np.cos(dd)*np.cos(rr), np.cos(dd)*np.sin(rr), np.sin(dd)])
    assert equal(rebuilt, xyz[ix])
    report['saved_core_chain'] = dict(cores=records, rows=len(ix), all_checks_exact=True)

    cache = np.load(OUT.parent / 'raw_environment_checks_20261014/raw_environment_catalogue.npz')
    choose = np.linspace(0, len(cache['TARGETID']) - 1, 256, dtype=int)
    rows = cache['catalogue_rows'][choose]
    product = read(ROOT / 'docs/figures/p12a_loa_cigale_20260929/PRODUCT.json')
    full = fitsio.read(product['path'], rows=rows); original = fitsio.read(vac, rows=rows)
    for name in original.dtype.names:
        assert equal(original[name], full[name]), name
    for name in ['TARGETID', 'RA', 'DEC', 'Z', 'P_WEB', 'EIGENVALUE_MEAN', 'CIGALE_MATCHED',
                 'MASS_CG_15', 'SFR_CG_15', 'AGE_CG_15', 'AV_CG_15', 'FAGN_CG_15']:
        cache_name = 'z' if name == 'Z' else name
        assert equal(cache[cache_name][choose], full[name]), name
    good = full['CIGALE_MATCHED']; selected = full[good]
    source = '/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits'
    c = fitsio.read(source, rows=selected['CIGALE_SOURCE_ROW'],
                    columns=['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'MASS_CG_15', 'SFR_CG_15', 'AGE_CG_15', 'AV_CG_15', 'FAGN_CG_15'])
    for name in ['TARGETID', 'MASS_CG_15', 'SFR_CG_15', 'AGE_CG_15', 'AV_CG_15', 'FAGN_CG_15']:
        assert equal(c[name], selected[name]), name
    dz = abs(c['Z'] - selected['Z']) / (1 + selected['Z'])
    sep = SkyCoord(c['TARGET_RA'], c['TARGET_DEC'], unit='deg').separation(SkyCoord(selected['RA'], selected['DEC'], unit='deg')).arcsec
    assert np.all(dz < .001) and np.all(sep < 1)
    report['plot_to_sed_chain'] = dict(sampled_rows=len(rows), cigale_rows=len(c), all_checks_exact=True,
                                      max_separation_arcsec=float(sep.max()), max_normalized_dz=float(dz.max()))
    report['script_sha256'] = sha(Path(__file__))
    (OUT / 'ordering_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
