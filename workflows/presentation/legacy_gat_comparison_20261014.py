#!/usr/bin/env python3
"""Saved GAT VAC audit and same-CIGALE comparison; never runs inference.

audit: bounded saved-catalogue inspection and legacy colour plots.
cigale: scan the 20 GB CIGALE source on an approved compute allocation,
or under an explicit bounded-login exception from the user.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import fitsio
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from shared.plot_style import apply_style, finalize_axes, CLASS_ORDER, COSMIC_WEB_COLORS
from desi_density_properties_20261014 import estimate

OLD = ROOT / 'DESI_BGS_PRERELEASE_VAC.pkl'
PRODUCT = ROOT / 'docs/figures/p12a_loa_cigale_20260929/PRODUCT.json'
SOURCE = Path('/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits')
CURRENT = Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010/raw_environment_checks_20261014/raw_environment_catalogue.npz')
NAMES = [x.title() for x in CLASS_ORDER]
COLORS = [COSMIC_WEB_COLORS[x] for x in CLASS_ORDER]
PROPS = ['MASS_CG_15', 'SFR_CG_15', 'AGE_CG_15', 'AV_CG_15', 'FAGN_CG_15']


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def save(fig, out, name):
    for ext in ['png', 'pdf']:
        fig.savefig(out / (name + '.' + ext), dpi=140, bbox_inches='tight')
    plt.close(fig)


def probabilities(d):
    p = d[['GAT_' + n.upper() + '_PROB' for n in NAMES]].to_numpy(float)
    assert np.isfinite(p).all() and (p >= 0).all() and (p <= 1).all()
    assert np.allclose(p.sum(1), 1, atol=1e-6)
    assert np.array_equal(p.argmax(1), d.GAT_ENV)
    return p


def colour_figure(mass, colour, p, out, name, title, colour_label='Observed g − r [mag]'):
    good = np.isfinite(mass) & np.isfinite(colour)
    mass, colour, p = mass[good], colour[good], p[good]
    hard = p.argmax(1)
    fig, ax = plt.subplots(2, 3, figsize=(16, 9), layout='constrained')
    fig.suptitle(title, fontsize=19)
    bins = np.linspace(0, 1.5, 101)
    for j in range(4):
        x = colour[hard == j]
        ax[0, 0].hist(x, bins=bins, histtype='step', color=COLORS[j], label=NAMES[j])
        ax[1, 0].hist(x, bins=bins, histtype='step', density=True, color=COLORS[j], label=NAMES[j])
    ax[0, 0].set_yscale('log')
    finalize_axes(ax[0, 0], 'Counts by label', colour_label, 'Galaxies')
    finalize_axes(ax[1, 0], 'Each class normalized separately', colour_label, 'Probability density')
    xb = np.linspace(9, 12, 61); yb = np.linspace(0, 1.5, 61)
    hs = [np.histogram2d(mass[hard == j], colour[hard == j], bins=[xb, yb])[0] for j in range(4)]
    # Fractions of all valid rows in each class, with one shared intensity scale.
    hs = [h / max(1, (hard == j).sum()) for j, h in enumerate(hs)]
    norm = LogNorm(vmin=1e-5, vmax=max(h.max() for h in hs))
    for j, a in enumerate(ax[:, 1:].flat):
        im = a.pcolormesh(xb, yb, np.ma.masked_less_equal(hs[j].T, 0), cmap='magma', norm=norm)
        finalize_axes(a, f'{NAMES[j]} · {(hard == j).sum():,}', 'log₁₀ stellar mass (catalogue units)', colour_label, legend=False)
    fig.colorbar(im, ax=ax[:, 1:], label='Fraction of class per bin', shrink=.7)
    save(fig, out, name)


def audit(out):
    d = pd.read_pickle(OLD); p = probabilities(d)
    assert d.TARGETID.dtype == np.dtype('int64')
    product = json.loads(PRODUCT.read_text())
    # Reading two columns only; preserve integer IDs throughout the exact join.
    v = fitsio.read(product['path'], columns=['TARGETID', 'Z'])
    common = np.intersect1d(d.TARGETID.to_numpy(), v['TARGETID'])
    overlaps = []
    for a in v[np.isin(v['TARGETID'], common)]:
        for _, b in d[d.TARGETID == a['TARGETID']].iterrows():
            overlaps.append({'TARGETID': int(a['TARGETID']), 'old_z': float(b.Z),
                             'new_z': float(a['Z']), 'dz_over_1plusz': float(abs(b.Z-a['Z'])/(1+a['Z']))})
    good = (d.FLUX_G > 0) & (d.FLUX_R > 0)
    colour = np.full(len(d), np.nan)
    colour[good] = -2.5 * np.log10(d.FLUX_G[good] / d.FLUX_R[good])
    r = {'legacy_path': str(OLD), 'legacy_sha256': sha(OLD), 'legacy_bytes': OLD.stat().st_size,
         'legacy_rows': len(d), 'unique_targetids': int(d.TARGETID.nunique()),
         'repeated_targetid_rows': int(d.TARGETID.duplicated(False).sum()),
         'repeated_targetids_with_different_labels': int((d.groupby('TARGETID').GAT_ENV.nunique()>1).sum()),
         'legacy_z_range': [float(d.Z.min()), float(d.Z.max())],
         'legacy_logmass_range': [float(d.LOGMSTAR.min()), float(d.LOGMSTAR.max())],
         'class_counts': dict(zip(NAMES, np.bincount(d.GAT_ENV, minlength=4).tolist())),
         'probability_checks_passed': True, 'median_max_probability': float(np.median(p.max(1))),
         'p12_product_receipt': product, 'p12_hash_freshly_computed': False,
         'p12_rows': len(v), 'p12_z_range': [float(v['Z'].min()), float(v['Z'].max())],
         'overlap_targetids': common.tolist(), 'overlap_count': len(common),
         'overlap_redshift_check': overlaps,
         'redshift_consistent_overlap_count': sum(x['dz_over_1plusz'] < .001 for x in overlaps),
         'colour_valid_rows': int(good.sum()),
         'property_origin': 'Legacy LOGMSTAR and fluxes are FastSpecFit catalogue columns; no CIGALE columns in pickle.',
         'colour_definition': '-2.5 log10(FLUX_G/FLUX_R); observed bands, no extinction or redshift correction',
         'exact_slide_provenance': 'Unconfirmed: saved pickle timestamp is January 2026, supplied slide is September 2025.',
         'script_sha256': sha(__file__)}
    for j, n in enumerate(NAMES):
        k = good & (d.GAT_ENV == j)
        r.setdefault('legacy_raw_means', {})[n] = {'mass': float(d.LOGMSTAR[k].mean()), 'colour': float(colour[k].mean())}
    dump(out / 'audit.json', r)
    colour_figure(d.LOGMSTAR.to_numpy(), colour, p, out, 'legacy_colour_mass',
                  'Saved GAT catalogue · 0.01 < z < 0.06 · FastSpecFit mass, observed colour')
    rest = (d.ABSMAG01_SDSS_G-d.ABSMAG01_SDSS_R).to_numpy()
    colour_figure(d.LOGMSTAR.to_numpy(), rest, p, out, 'legacy_rest_colour_mass',
                  'Saved GAT catalogue · FastSpecFit mass and rest-frame colour',
                  'SDSS g − r at reference z = 0.1 [mag]')
    k = np.isfinite(rest) & ~d.TARGETID.duplicated(False).to_numpy() & (d.LOGMSTAR >= 9) & (d.LOGMSTAR < 12)
    y = rest[k]; m = d.LOGMSTAR.to_numpy()[k]; z = d.Z.to_numpy()[k]
    w = np.eye(4)[p[k].argmax(1)]
    cells = np.floor((m-9)/.1).astype(int)*10 + np.floor(z/.01).astype(int)
    _, cells = np.unique(cells, return_inverse=True)
    n = np.stack([np.bincount(cells, weights=w[:, j]) for j in range(4)])
    s = np.stack([np.bincount(cells, weights=w[:, j]*y) for j in range(4)])
    keep = (n >= 20).all(0); ref = n[:, keep].sum(0); ref /= ref.sum()
    raw = w.T@y/w.sum(0); controlled = (s[:, keep]/n[:, keep]*ref).sum(1)
    dump(out/'legacy_colour_control.json', {'definition': 'ABSMAG01_SDSS_G - ABSMAG01_SDSS_R; FastSpecFit log mass as saved',
         'unique_only_rows': int(k.sum()), 'controlled_rows': int(keep[cells].sum()), 'cells': int(keep.sum()),
         'raw_means': raw.tolist(), 'controlled_means': controlled.tolist(),
         'raw_cluster_minus_void': float(raw[-1]-raw[0]), 'controlled_cluster_minus_void': float(controlled[-1]-controlled[0]),
         'control': '.1 dex saved mass and .01 z, >=20 per hard class, pooled fixed common cell mixture',
         'uncertainty': 'Not estimated in this bounded diagnostic'})
    src = Path('/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/loa-combined-lowz.fits')
    c = fitsio.read(src, columns=['TARGETID', 'ZWARN', 'DELTACHI2', 'LOGMSTAR', 'BGS_TARGET'])
    selected = (c['ZWARN']==0)&(c['DELTACHI2']>=25)&(c['LOGMSTAR']>=9)&(c['BGS_TARGET']!=0)
    dump(out/'selection_check.json', {'source': str(src), 'source_rows': len(c), 'simple_cut_rows': int(selected.sum()),
         'legacy_rows': len(d), 'same_targetid_multiset': bool(np.array_equal(np.sort(c['TARGETID'][selected]), np.sort(d.TARGETID.to_numpy()))),
         'cuts': 'ZWARN==0, DELTACHI2>=25, LOGMSTAR>=9, BGS_TARGET!=0',
         'extra_mass_completeness_filter_visible_in_current_loader': False})
    print(json.dumps({k: r[k] for k in ['legacy_rows', 'legacy_z_range', 'p12_z_range', 'overlap_count', 'class_counts', 'median_max_probability', 'legacy_raw_means']}, indent=2))


def join_cigale(out):
    d = pd.read_pickle(OLD)
    # Do not choose between repeated spectra with different saved predictions.
    duplicate_rows = int(d.TARGETID.duplicated(False).sum())
    d = d[~d.TARGETID.duplicated(False)].sort_values('TARGETID').reset_index(drop=True)
    ids = d.TARGETID.to_numpy(); probabilities(d)
    cols = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'SPECTYPE', 'FLUX_G', 'FLUX_R'] + PROPS
    parts = []; source_rows = []; id_matches = 0
    with fitsio.FITS(SOURCE) as f:
        for lo in range(0, f[1].get_nrows(), 500000):
            c = f[1].read(rows=np.arange(lo, min(lo + 500000, f[1].get_nrows())), columns=cols)
            at = np.searchsorted(ids, c['TARGETID']); at = np.minimum(at, len(ids)-1)
            keep = ids[at] == c['TARGETID']; id_matches += int(keep.sum())
            if keep.any():
                parts.append(c[keep]); source_rows.append(lo + np.flatnonzero(keep))
            if lo % 5000000 == 0: print('CIGALE rows scanned:', lo, flush=True)
    c = np.concatenate(parts); sr = np.concatenate(source_rows); at = np.searchsorted(ids, c['TARGETID'])
    dz = np.abs(c['Z']-d.Z.to_numpy()[at]) / (1+d.Z.to_numpy()[at])
    dra = (c['TARGET_RA']-d.RA.to_numpy()[at]+180) % 360-180
    sep = 3600*np.hypot(dra*np.cos(np.deg2rad(d.DEC.to_numpy()[at])), c['TARGET_DEC']-d.DEC.to_numpy()[at])
    ok = np.isfinite(dz) & np.isfinite(sep) & (dz < .001) & (sep < 1)
    rejected = int((~ok).sum()); c, sr, at = c[ok], sr[ok], at[ok]
    unique, counts = np.unique(at, return_counts=True)
    keep = np.isin(at, unique[counts == 1]); ambiguous = int((counts > 1).sum())
    c, sr, at = c[keep], sr[keep], at[keep]
    old = {k: c[k] for k in ['TARGETID', 'SPECTYPE', 'FLUX_G', 'FLUX_R'] + PROPS}
    old.update(z=d.Z.to_numpy()[at], RA=d.RA.to_numpy()[at], DEC=d.DEC.to_numpy()[at],
               P_WEB=probabilities(d)[at], source_row=sr)
    np.savez_compressed(out/'legacy_cigale.npz', **old)
    dump(out/'join.json', {'source': str(SOURCE), 'source_bytes': SOURCE.stat().st_size,
         'source_hash_from_existing_receipt': json.loads((ROOT/'docs/figures/p12a_loa_cigale_20260929/JOIN.json').read_text())['source_sha256'],
         'source_hash_freshly_computed': False, 'legacy_sha256': sha(OLD), 'id_matched_source_rows': id_matches,
         'coordinate_redshift_rejected_rows': rejected, 'ambiguous_ids_excluded': ambiguous,
         'accepted_unique_ids': len(c), 'output_sha256': sha(out/'legacy_cigale.npz'),
         'legacy_repeated_targetid_rows_excluded': duplicate_rows,
         'rule': 'Exact integer TARGETID, separation <1 arcsec, abs(dz)/(1+z)<.001, only one consistent observation.'})


def compare(out):
    import healpy as hp
    datasets = {'GAT': dict(np.load(out/'legacy_cigale.npz')), 'P12-A': dict(np.load(CURRENT))}
    provenance = json.loads(CURRENT.with_name('measurement_provenance.json').read_text())
    assert sha(CURRENT) == provenance['catalogue_sha256']
    stats = {}; data = {}; apply_style()
    for model, a in datasets.items():
        lm = np.log10(np.maximum(a['MASS_CG_15'], 1e-30)); ls = np.log10(np.maximum(a['SFR_CG_15'], 1e-30))
        good = np.ones(len(lm), bool)
        if model == 'P12-A': good &= a['CIGALE_MATCHED']; spec = a['CIGALE_SPECTYPE']
        else: spec = a['SPECTYPE']
        good &= np.char.strip(spec.astype('U')) == 'GALAXY'
        for prop in PROPS[:3]: good &= np.isfinite(a[prop]) & (a[prop] > 0)
        good &= np.isfinite(a['AV_CG_15']) & (a['AV_CG_15'] >= 0)
        good &= np.isfinite(a['FAGN_CG_15']) & (a['FAGN_CG_15'] >= 0) & (a['FAGN_CG_15'] <= 1)
        good &= (lm >= 9) & (lm < 12)
        y = {'mass': lm, 'sfr': ls, 'ssfr': ls-lm, 'quenched': (ls-lm < -11).astype(float),
             'age': np.log10(np.maximum(a['AGE_CG_15'],1e-30)/1000), 'av': a['AV_CG_15'], 'fagn': a['FAGN_CG_15']}
        p = a['P_WEB'][good].astype(float); assert np.allclose(p.sum(1),1,atol=1e-6)
        weights = {'probability': p, 'label': np.eye(4)[p.argmax(1)]}
        reg = hp.ang2pix(8, a['RA'][good], a['DEC'][good], lonlat=True)
        dz = .01; zb = np.floor(a['z'][good]/dz).astype(int)
        mb = np.floor((lm[good]-9)/.1).astype(int)
        data[model] = {'mass': lm[good], 'ssfr': y['ssfr'][good], 'quenched': y['quenched'][good],
                       'p': p, 'weights': weights, 'region': reg, 'mb': mb, 'zb': zb}
        stats[model] = {'valid_rows': int(good.sum()), 'median_max_probability': float(np.median(p.max(1))), 'properties': {}}
        for name, value in y.items():
            yy = value[good]; cells = zb if name == 'mass' else mb*100+zb
            result, keep = estimate(yy, weights, cells, reg, np.zeros(len(yy)), nboot=64)
            result.pop('paired_density_minus_soft_contrast')
            result['raw_means'] = {n: (w.T@yy/w.sum(0)).tolist() for n,w in weights.items()}
            # Same retained rows for the raw-v-controlled comparison.
            result['raw_same_rows'] = {n: (w[keep].T@yy[keep]/w[keep].sum(0)).tolist() for n,w in weights.items()}
            stats[model]['properties'][name] = result
        if model == 'GAT':
            color = np.full(len(lm), np.nan); k = (a['FLUX_G']>0)&(a['FLUX_R']>0)
            color[k] = -2.5*np.log10(a['FLUX_G'][k]/a['FLUX_R'][k])
            colour_figure(lm[good], color[good], a['P_WEB'][good], out, 'legacy_cigale_colour_mass',
                          'GAT · current CIGALE mass and observed colour · 0.01 < z < 0.06')
    dump(out/'comparison.json', {'samples': stats, 'current_cache': str(CURRENT), 'current_sha256': sha(CURRENT),
         'script_sha256': sha(__file__), 'estimator_sha256': sha(Path(__file__).with_name('desi_density_properties_20261014.py')),
         'legacy_join_sha256': sha(out/'legacy_cigale.npz'),
         'limits': 'Disjoint redshift ranges and sky samples. Each sample separately balanced in mass/z; not a same-galaxy or common-redshift model test. No environment truth or causal claim.',
         'controls': '.1 dex mass and .01 z; at least 20 effective rows per class for both label/probability methods. Fixed common mixture within each sample. Mass controls z only.',
         'uncertainty': '16-84% intervals from 64 nside8 sky-block resamples; excludes property/model/selection uncertainty.',
         'property_notes': 'CG15; log sSFR<-11; fitted AGE, not established mass-weighted age; FAGN is fitted contribution, not host incidence.'})
    labels = [r'log₁₀ stellar mass [$M_\odot$]', r'log₁₀ SFR [$M_\odot$/yr]', 'log₁₀ sSFR [1/yr]', 'Quenched fraction', 'log₁₀ fitted age [Gyr]', r'Dust $A_V$ [mag]', 'AGN contribution']
    for mode in ['raw_same_rows', 'controlled']:
        fig, axs = plt.subplots(2,4,figsize=(18,9),layout='constrained')
        fig.suptitle('Current CIGALE properties · ' + ('raw means on retained rows' if mode=='raw_same_rows' else 'mass and redshift balanced within each sample'),fontsize=18)
        for ax, key, label in zip(axs.flat, stats['GAT']['properties'], labels):
            for model,col in [('GAT','#F5C144'),('P12-A','#4E84F7')]:
                r = stats[model]['properties'][key]
                for method,ls,marker in [('probability','-','o'),('label','--','s')]:
                    mu = r['raw_same_rows'][method] if mode=='raw_same_rows' else r['methods'][method]['mean']
                    ax.plot(range(4),mu,color=col,ls=ls,marker=marker,label=model+' '+method)
                    if mode=='controlled':
                        ci = r['methods'][method]['interval']; ax.vlines(range(4),ci[0],ci[1],color=col,lw=1)
            ax.set_xticks(range(4),['V','W','F','C'])
            finalize_axes(ax,label,'Environment','Mean' if key!='quenched' else 'Fraction',legend=False)
            if key=='quenched': ax.set_ylim(0,1)
        axs.flat[-1].axis('off'); axs.flat[-1].legend(*axs.flat[0].get_legend_handles_labels(),loc='upper left')
        axs.flat[-1].text(0,.48,'GAT: z = 0.01–0.06\nP12-A: z = 0.20–0.30 wedge\n\nDifferent galaxies and epochs.\nThis does not isolate the model.\n\nV/W/F/C: void/wall/filament/cluster',va='top',fontsize=11)
        save(fig,out,'cigale_properties_'+mode)
    distribution_figure(data, out)
    common_mass_check(data, out)
    contrast_figure(stats, out)


def distribution_figure(data, out):
    xb = np.linspace(9,12,61); yb = np.linspace(-14,-8,61)
    hist = {}
    for model,a in data.items():
        hist[model] = [np.histogram2d(a['mass'],a['ssfr'],bins=[xb,yb],weights=a['p'][:,j])[0]/a['p'][:,j].sum() for j in range(4)]
    norm = LogNorm(1e-5, max(h.max() for hs in hist.values() for h in hs))
    fig,axs = plt.subplots(2,4,figsize=(17,8),layout='constrained',sharex=True,sharey=True)
    fig.suptitle('Same CIGALE measurements · different galaxy samples',fontsize=20)
    for i,(model,hs) in enumerate(hist.items()):
        for j,h in enumerate(hs):
            ax = axs[i,j]; im = ax.pcolormesh(xb,yb,np.ma.masked_less_equal(h.T,0),cmap='magma',norm=norm)
            ax.axhline(-11,color='white',lw=.7,ls=':')
            finalize_axes(ax,model+' · '+NAMES[j], r'log₁₀ stellar mass [$M_\odot$]', 'log₁₀ sSFR [1/yr]',legend=False)
    fig.colorbar(im,ax=axs,label='Probability-weighted fraction of class per bin',shrink=.7)
    save(fig,out,'cigale_mass_ssfr_distributions')


def common_mass_check(data, out):
    """Same mass mixture across samples; each keeps its own z distribution.

    All label/soft groups share support. Within each retained mass bin, use
    each sample's pooled eligible redshift mixture. Across mass bins use the
    smaller pooled count of the two samples, normalized to sum to one.
    """
    tables = {}; eligible_mass = []
    for model,a in data.items():
        ids,c = np.unique(a['mb']*100+a['zb'],return_inverse=True)
        _,reg = np.unique(a['region'],return_inverse=True); nr=reg.max()+1; nc=len(ids)
        flat=c*nr+reg; weights=np.concatenate(list(a['weights'].values()),axis=1)
        count=np.stack([np.bincount(flat,weights=weights[:,j],minlength=nc*nr).reshape(nc,nr) for j in range(8)])
        square=np.stack([np.bincount(c,weights=weights[:,j]**2,minlength=nc) for j in range(8)])
        good=((count.sum(-1)**2/np.maximum(square,1e-30))>=20).all(0)
        sums={key:np.stack([np.bincount(flat,weights=weights[:,j]*a[key],minlength=nc*nr).reshape(nc,nr) for j in range(8)]) for key in ['ssfr','quenched']}
        pooled=np.bincount(c,minlength=nc); mass=ids//100
        tables[model]={'count':count,'sums':sums,'good':good,'pooled':pooled,'mass':mass,'rowcells':c,'nr':nr}
        eligible_mass.append(np.unique(mass[good]))
    common=np.intersect1d(*eligible_mass)
    assert len(common)>0
    totals=[]
    for t in tables.values():
        t['good'] &= np.isin(t['mass'],common)
        totals.append(np.array([t['pooled'][t['good']&(t['mass']==b)].sum() for b in common]))
    mass_ref=np.minimum(*totals).astype(float); mass_ref/=mass_ref.sum()
    results={}
    for model,t in tables.items():
        k=t['good']; m=t['mass'][k]; n=t['pooled'][k]; ref=np.zeros(len(m))
        for b,prob in zip(common,mass_ref):
            sel=m==b; ref[sel]=prob*n[sel]/n[sel].sum()
        assert np.isclose(ref.sum(),1)
        count=t['count'][:,k]; sums={key:v[:,k] for key,v in t['sums'].items()}
        mean={key:(v.sum(-1)/count.sum(-1)*ref).sum(-1) for key,v in sums.items()}
        rng=np.random.default_rng(20261014); boot={key:[] for key in sums}
        for _ in range(64):
            mult=np.bincount(rng.integers(t['nr'],size=t['nr']),minlength=t['nr']); den=count@mult
            if (den>0).all():
                for key,v in sums.items():boot[key].append(((v@mult)/den*ref).sum(-1))
        assert min(map(len,boot.values()))>=52
        result={'rows':int(k[t['rowcells']].sum()),'cells':int(k.sum()),'properties':{}}
        for key,mu in mean.items():
            bs=np.array(boot[key]); result['properties'][key]={}
            for j,method in enumerate(['probability','label']):
                offset=4*j; result['properties'][key][method]={'mean':mu[offset:offset+4].tolist(),
                  'contrast':float(mu[offset+3]-mu[offset]),
                  'interval':np.quantile(bs[:,offset+3]-bs[:,offset],[.16,.84]).tolist()}
        results[model]=result
    dump(out/'common_mass_check.json', {'mass_bin_lower_edges':(9+.1*common).tolist(),'mass_reference_weights':mass_ref.tolist(),
         'results':results,'definition':common_mass_check.__doc__,
         'limits':'Matching mass cannot match the disjoint redshift ranges or change the graph tracer populations. Sky coverage, SED and model-transfer differences remain.'})


def contrast_figure(stats, out):
    matched=json.loads((out/'common_mass_check.json').read_text())['results']
    fig,axs=plt.subplots(1,2,figsize=(14,7.5))
    fig.subplots_adjust(left=.075,right=.98,bottom=.27,top=.83,wspace=.27)
    fig.suptitle('Matching stellar mass and label treatment narrows the gap',fontsize=18,y=.96)
    for ax,key,factor,title in zip(axs,['quenched','ssfr'],[100,-1],
                                  ['Cluster − void quenched fraction','Void − cluster mean log sSFR']):
        for model,color in [('GAT','#F5C144'),('P12-A','#4E84F7')]:
            first=stats[model]['properties'][key]['methods']['probability']
            second=matched[model]['properties'][key]['probability']
            third=matched[model]['properties'][key]['label']
            points=[first['contrast'],second['contrast'],third['contrast']]
            intervals=[first['contrast_interval'],second['interval'],third['interval']]
            points=np.array(points)*factor; intervals=np.sort(np.array(intervals)*factor,axis=1)
            offset=-.045 if model=='GAT' else .045; x=np.arange(3)+offset
            ax.errorbar(x,points,yerr=[points-intervals[:,0],intervals[:,1]-points],fmt='o',color=color,capsize=4,label=model)
            for xx,yy in zip(x,points):ax.annotate(f'{yy:.1f}' if key=='quenched' else f'{yy:.3f}',(xx,yy),xytext=(7,6),textcoords='offset points',fontsize=11,color=color)
        ax.set_xticks(range(3),['Each sample’s\nmass mixture','Shared\nmass mixture','Shared mass mixture\n+ hard labels'])
        ax.set_xlim(-.3,2.45); ax.set_ylim(bottom=0)
        finalize_axes(ax,title,'','Percentage points' if key=='quenched' else 'dex',legend=True)
        ax.grid(axis='y',alpha=.15)
    fig.text(.5,.045,'All stages balance redshift within each sample. Shared mass: 10.2 ≤ log₁₀ M < 11.3.\nGAT z=0.01–0.06; P12-A z=0.20–0.30. Different samples remain; this is not a model-accuracy test.\nIntervals: 16–84% sky-block resampling; stages change retained rows and weighting.',ha='center',fontsize=11,linespacing=1.6)
    save(fig,out,'selection_and_weighting_check')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--stage',choices=['audit','cigale','compare'],required=True); args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True); apply_style()
    if args.stage=='audit': audit(args.out)
    if args.stage=='cigale': join_cigale(args.out); compare(args.out)
    if args.stage=='compare': compare(args.out)
