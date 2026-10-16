"""Paired CG5/CG15 check on the existing P12-A wedge; no new inference."""
import json
import sys
from pathlib import Path

import fitsio
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from shared.plot_style import apply_style, finalize_axes
from desi_density_properties_20261014 import estimate, sha, dump

OUT = Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010/p12a_cg5_check_20261014')
CACHE = OUT.parent / 'raw_environment_checks_20261014/raw_environment_catalogue.npz'
LABELS = {'mass': r'log$_{10}$ stellar mass [$M_\odot$]',
          'sfr': r'log$_{10}$ SFR [$M_\odot$/yr]', 'ssfr': r'log$_{10}$ sSFR [1/yr]',
          'quenched': 'Quenched fraction', 'age': r'log$_{10}$ fitted age [Gyr]',
          'av': r'Dust $A_V$ [mag]', 'fagn': 'Fitted AGN contribution'}


def properties(d, fit):
    m, s, age, av, agn = [d[f'{p}_CG_{fit}'] for p in ['MASS', 'SFR', 'AGE', 'AV', 'FAGN']]
    good = np.isfinite([m, s, age, av, agn]).all(0)
    good &= (m > 0) & (s > 0) & (age > 0) & (av >= 0) & (agn >= 0) & (agn <= 1)
    lm = np.log10(np.maximum(m, 1e-30)); ls = np.log10(np.maximum(s, 1e-30))
    good &= (lm >= 9) & (lm < 12)
    return good, dict(mass=lm, sfr=ls, ssfr=ls-lm, quenched=(ls-lm < -11).astype(float),
                      age=np.log10(np.maximum(age, 1e-30)/1000), av=av, fagn=agn)


def main():
    OUT.mkdir(exist_ok=True)
    receipt = json.loads((ROOT/'docs/figures/p12a_loa_cigale_20260929/PRODUCT.json').read_text())
    a = np.load(CACHE); provenance = json.loads(CACHE.with_name('measurement_provenance.json').read_text())
    assert sha(CACHE) == provenance['catalogue_sha256']
    cols = ['TARGETID', 'Z', 'SUPPORTED', 'P_WEB', 'CIGALE_MATCHED', 'CIGALE_SPECTYPE', 'CIGALE_SOURCE_ROW']
    cols += [f'{p}_CG_{f}' for f in [15, 5] for p in ['MASS', 'SFR', 'AGE', 'AV', 'FAGN']]
    d = fitsio.read(receipt['path'], rows=a['catalogue_rows'], columns=cols)
    assert np.array_equal(d['TARGETID'], a['TARGETID']) and np.array_equal(d['Z'], a['z'])
    assert np.array_equal(d['P_WEB'], a['P_WEB'], equal_nan=True)
    for p in ['MASS', 'SFR', 'AGE', 'AV', 'FAGN']:
        assert np.array_equal(d[f'{p}_CG_15'], a[f'{p}_CG_15'], equal_nan=True)
    base = d['SUPPORTED'] & d['CIGALE_MATCHED'] & (np.char.strip(d['CIGALE_SPECTYPE'].astype('U')) == 'GALAXY')
    g15, y15 = properties(d, 15); g5, y5 = properties(d, 5); keep = base & g15 & g5
    p = d['P_WEB'][keep].astype(float); weights = {'probability': p, 'label': np.eye(4)[p.argmax(1)]}
    zbin = np.floor(d['Z'][keep]/.01).astype(int)
    cells15 = zbin + 100*np.floor((y15['mass'][keep]-9)/.1).astype(int)
    cells5 = zbin + 100*np.floor((y5['mass'][keep]-9)/.1).astype(int)
    region = a['region'][keep]; zeros = np.zeros(keep.sum()); reports = {}
    for name in LABELS:
        cells = zbin if name == 'mass' else cells15
        result = {}
        for fit, yy in [('CG15', y15[name][keep]), ('CG5', y5[name][keep]), ('CG5_minus_CG15', y5[name][keep]-y15[name][keep])]:
            r, support = estimate(yy, weights, cells, region, zeros, nboot=64)
            r.pop('paired_density_minus_soft_contrast')
            r['raw_same_rows'] = {k: (w[support].T @ yy[support]/w[support].sum(0)).tolist() for k, w in weights.items()}
            result[fit] = r
        assert result['CG15']['common_cell_ids'] == result['CG5']['common_cell_ids']
        assert result['CG15']['reference_weights'] == result['CG5']['reference_weights']
        for method in weights:
            expected = np.asarray(result['CG5']['methods'][method]['mean']) - result['CG15']['methods'][method]['mean']
            assert np.allclose(expected, result['CG5_minus_CG15']['methods'][method]['mean'], atol=1e-12)
        reports[name] = result
        print(name, {fit: round(result[fit]['methods']['probability']['contrast'], 5) for fit in result}, flush=True)
    sensitivity = {}
    for name in ['ssfr', 'quenched']:
        r, _ = estimate(y5[name][keep], weights, cells5, region, zeros, nboot=64)
        r.pop('paired_density_minus_soft_contrast'); sensitivity[name] = r
    # Fresh direct source reads verify CG5 was carried through the join correctly.
    ix = np.flatnonzero(keep)[np.linspace(0, keep.sum()-1, 64, dtype=int)]
    source = '/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits'
    checkcols = ['TARGETID'] + [f'{p}_CG_{f}' for f in [15, 5] for p in ['MASS', 'SFR', 'AGE', 'AV', 'FAGN']]
    c = fitsio.read(source, rows=d['CIGALE_SOURCE_ROW'][ix], columns=checkcols)
    for name in checkcols:
        assert np.array_equal(c[name], d[name][ix], equal_nan=True), name
    out = dict(rows={'wedge':len(d), 'CG15_valid':int((base & g15).sum()), 'CG5_valid':int((base & g5).sum()),
                     'joint_valid':int(keep.sum())}, properties=reports, CG5_mass_control_sensitivity=sensitivity,
               method='Identical galaxies/probabilities; .1dex CG15 mass/.01z cells, pooled common weights, >=20 effective rows in every hard/soft class. Mass uses z only. CG5 own-mass sensitivity separately reported.',
               uncertainty='64 paired nside8 sky-block resamples; 16-84% intervals. CG5-minus-CG15 intervals use per-galaxy differences with the same resampling and weights.',
               limits='Existing supported z.20-.30 wedge. Same CIGALE source and redshift, not independent SED truth; photometric/model errors and selection uncertainty excluded.',
               cache_sha256=sha(CACHE), script_sha256=sha(__file__), source_receipt=receipt,
               direct_source_rows_verified=len(ix), no_new_inference=True)
    dump(OUT/'RESULTS.json', out)
    apply_style()
    fig, axes = plt.subplots(2, 4, figsize=(17, 8.5), layout='constrained')
    fig.suptitle('P12-A environments: CG5 versus CG15 on the same galaxies', fontsize=18)
    for ax, (name, label) in zip(axes.flat, LABELS.items()):
        for fit, color in [('CG15', '#F5C144'), ('CG5', '#4E84F7')]:
            for method, ls, marker in [('probability', '-', 'o'), ('label', '--', 's')]:
                r = reports[name][fit]['methods'][method]
                ax.plot(range(4), r['mean'], color=color, ls=ls, marker=marker, label=fit+' '+method)
                ax.vlines(range(4), *r['interval'], color=color, lw=.8)
        ax.set_xticks(range(4), ['Void', 'Wall', 'Fil.', 'Cl.'])
        finalize_axes(ax, label, 'Environment', 'Mean', legend=False)
        if name == 'quenched': ax.set_ylim(0, 1)
    axes.flat[-1].axis('off')
    axes.flat[-1].legend(*axes.flat[0].get_legend_handles_labels(), loc='upper left')
    axes.flat[-1].text(0, .52, f"{keep.sum():,} jointly valid galaxies\nz = 0.20–0.30 wedge\n\nSame CG15 mass/redshift controls.\nMass panel controls redshift only.\n\nBars: 16–84% spatial resampling.\nFits share observations and assumptions.", va='top', fontsize=11)
    for ext in ['png', 'pdf']: fig.savefig(OUT/f'cg5_vs_cg15.{ext}', dpi=140, bbox_inches='tight')
    plt.close(fig)
    print('COMPLETE', OUT, flush=True)


if __name__ == '__main__':
    main()
