"""Sequential fixed-prediction selection/confidence checks; CPU analysis only.

prepare scans the registered 8 GB Loa targeting catalogue once, in 250k-row
blocks. analyse uses the resulting small cache. No training or model inference.
"""
import argparse
import json
from pathlib import Path

import fitsio
import healpy as hp
import numpy as np
import pandas as pd

from desi_density_properties_20261014 import estimate, sha, dump

ROOT = Path('/global/homes/d/dkololgi/GraphWeb_DESI')
PRESENTATION = Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010')
OLD = PRESENTATION/'legacy_gat_comparison_20261014'
CURRENT = PRESENTATION/'raw_environment_checks_20261014/raw_environment_catalogue.npz'
OUT = PRESENTATION/'legacy_selection_tests_20261014'


def sky(ra, dec):
    return (ra >= 120) & (ra < 160) & (dec >= 14.5) & (dec < 30.6)


def completeness(z, mass, r, rlimit, edges):
    """Within each z slice, faintest 20%, then 90th-percentile limiting mass."""
    bins = np.searchsorted(edges, z, side='right')-1
    result = []
    for i, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        k = (bins == i) & np.isfinite(mass) & np.isfinite(r) & (r < rlimit)
        if k.sum() < 100:
            result.append(dict(zlo=float(lo), zhi=float(hi), rows=int(k.sum()), usable=False)); continue
        # Rank magnitude relative to each object's actual targeting limit.
        offset = r-rlimit; faint = k & (offset >= np.quantile(offset[k], .8))
        limits = mass[faint] + .4*(r[faint]-rlimit[faint])
        result.append(dict(zlo=float(lo), zhi=float(hi), rows=int(k.sum()), faint_rows=int(faint.sum()),
                           logmass90=float(np.quantile(limits, .9)), logmass95=float(np.quantile(limits, .95)), usable=True))
    return result


def apply_limit(z, mass, curve, extra=0.):
    keep = np.zeros(len(z), dtype=bool)
    for b in curve:
        if b['usable']:
            keep |= (z >= b['zlo']) & (z < b['zhi']) & (mass >= b['logmass90']+extra)
    return keep


def prepare():
    OUT.mkdir(exist_ok=True)
    current = dict(np.load(CURRENT)); ids = current['TARGETID']; order = np.argsort(ids); sid = ids[order]
    ready = json.loads(Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1/INPUTS_READY.json').read_text())
    source = ready['source']['path']; r = np.full(len(ids), np.nan); lim = r.copy(); found = np.zeros(len(ids), int)
    photsys = np.full(len(ids), '', dtype='U1')
    with fitsio.FITS(source) as f:
        for start in range(0, f[1].get_nrows(), 250000):
            d = f[1].read(rows=np.arange(start, min(start+250000, f[1].get_nrows())),
                          columns=['TARGETID', 'RA', 'DEC', 'Z_not4clus', 'FLUX_R', 'MW_TRANSMISSION_R', 'PHOTSYS'])
            at = np.searchsorted(sid, d['TARGETID']); ok = (at < len(ids)) & (sid[np.minimum(at, len(ids)-1)] == d['TARGETID'])
            if ok.any():
                ii = order[at[ok]]; dd = d[ok]
                for name, values in [('RA', current['RA']), ('DEC', current['DEC']), ('Z_not4clus', current['z'])]:
                    assert np.array_equal(dd[name], values[ii]), name
                flux = dd['FLUX_R']/dd['MW_TRANSMISSION_R']
                assert np.all(flux > 0)
                r[ii] = 22.5-2.5*np.log10(flux); photsys[ii] = np.char.strip(dd['PHOTSYS'].astype('U'))
                lim[ii] = np.where(photsys[ii] == 'N', 19.54, 19.5); found[ii] += 1
    assert np.all(found == 1)
    np.savez_compressed(OUT/'p12a_targeting.npz', TARGETID=ids, r=r, rlimit=lim, photsys=photsys)
    dump(OUT/'PREPARATION.json', dict(source=ready['source'], full_source_hash_recomputed=False,
          cache_sha256=sha(CURRENT), targeting_sha256=sha(OUT/'p12a_targeting.npz'), rows=len(ids),
          exact_id_sky_redshift_match=True, r_definition='22.5-2.5log10(FLUX_R/MW_TRANSMISSION_R)',
          limits='Installed desitarget.cuts.isBGS_colors: BRIGHT South19.5, North19.54. Main BRIGHT comparison only.',
          photsys_counts={k:int((photsys==k).sum()) for k in np.unique(photsys)}))
    print('PREPARED', len(ids), flush=True)


def valid(a, current=False):
    good = (a['CIGALE_MATCHED'].copy() if current else np.ones(len(a['z']), bool))
    good &= np.char.strip(a['CIGALE_SPECTYPE' if current else 'SPECTYPE'].astype('U')) == 'GALAXY'
    for prop in ['MASS_CG_15', 'SFR_CG_15', 'AGE_CG_15']:
        good &= np.isfinite(a[prop]) & (a[prop] > 0)
    good &= np.isfinite(a['AV_CG_15']) & (a['AV_CG_15'] >= 0)
    good &= np.isfinite(a['FAGN_CG_15']) & (a['FAGN_CG_15'] >= 0) & (a['FAGN_CG_15'] <= 1)
    mass = np.log10(np.maximum(a['MASS_CG_15'], 1e-30)); ssfr = np.log10(np.maximum(a['SFR_CG_15'],1e-30))-mass
    good &= (mass >= 9) & (mass < 12)
    return good, mass, ssfr


def paired_selection(a, base, subset, mass, ssfr, label):
    """Same mass/z cells and reference mixture for parent and selected samples."""
    ii = np.flatnonzero(base); sub = subset[ii]
    p = a['P_WEB'][ii].astype(float); hard = np.eye(4)[p.argmax(1)]
    weights = {'all_soft':p, 'all_hard':hard, 'selected_soft':p*sub[:,None], 'selected_hard':hard*sub[:,None]}
    zbin = np.floor(a['z'][ii]/.01).astype(int); mb = np.floor((mass[ii]-9)/.1).astype(int); cells = mb*100+zbin
    reg = hp.ang2pix(8, a['RA'][ii], a['DEC'][ii], lonlat=True)
    result = dict(definition=label, parent_rows=len(ii), selected_rows=int(sub.sum()), properties={})
    try:
        for name, yy in [('ssfr', ssfr[ii]), ('quenched', (ssfr[ii] < -11).astype(float))]:
            r, keep = estimate(yy, weights, cells, reg, np.zeros(len(ii)), minimum=20, nboot=64)
            r.pop('paired_density_minus_soft_contrast')
            r['selected_retained_rows'] = int((keep & sub).sum())
            result['properties'][name] = r
        result['status'] = 'measured'
    except (ValueError, AssertionError) as e:
        result['status'] = 'insufficient_common_support'; result['reason'] = str(e)
    return result


def analyse():
    old = dict(np.load(OLD/'legacy_cigale.npz')); new = dict(np.load(CURRENT))
    receipt = json.loads((OUT/'PREPARATION.json').read_text())
    assert sha(CURRENT) == receipt['cache_sha256'] and sha(OUT/'p12a_targeting.npz') == receipt['targeting_sha256']
    phot = np.load(OUT/'p12a_targeting.npz'); assert np.array_equal(phot['TARGETID'],new['TARGETID'])
    saved = pd.read_pickle(ROOT/'DESI_BGS_PRERELEASE_VAC.pkl')
    saved = saved[~saved.TARGETID.duplicated(False)].set_index('TARGETID').loc[old['TARGETID']]
    assert np.array_equal(saved.Z.to_numpy(),old['z'])
    old_good, old_mass, old_ssfr = valid(old); new_good, new_mass, new_ssfr = valid(new,True)
    bright = (saved.BGS_TARGET.to_numpy() & 2) != 0
    main = np.char.strip(saved.SURVEY.to_numpy().astype('U')) == 'main'
    footprint = sky(old['RA'],old['DEC'])
    old_r = 22.5-2.5*np.log10(np.maximum(saved.FLUX_R.to_numpy(),1e-30))
    old_base = old_good & bright & main & footprint & (old_r < 19.5)
    new_base = new_good & (phot['r'] < phot['rlimit'])
    # Derive the old completeness curve BEFORE the saved VAC's mass cut.
    parent = fitsio.read('/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/loa-combined-lowz.fits',
                        columns=['TARGETID','SURVEY','RA','DEC','Z','ZWARN','DELTACHI2','BGS_TARGET','FLUX_R','LOGMSTAR'])
    _, inv, count = np.unique(parent['TARGETID'],return_inverse=True,return_counts=True)
    pk = (count[inv]==1) & (parent['ZWARN']==0) & (parent['DELTACHI2']>=25) & ((parent['BGS_TARGET'] & 2)!=0)
    pk &= (np.char.strip(parent['SURVEY'].astype('U'))=='main') & sky(parent['RA'],parent['DEC']) & (parent['FLUX_R']>0)
    pm = parent['LOGMSTAR'][pk]-2*np.log10(.704); pr = 22.5-2.5*np.log10(parent['FLUX_R'][pk])
    old_curve = completeness(parent['Z'][pk],pm,pr,np.full(len(pr),19.5),np.linspace(.01,.06000001,11))
    # No stellar-mass truncation in the current parent for estimating the limit.
    nk = new['CIGALE_MATCHED'] & np.isfinite(new_mass) & (new['MASS_CG_15']>0) & (phot['r']<phot['rlimit'])
    new_curve = completeness(new['z'][nk],new_mass[nk],phot['r'][nk],phot['rlimit'][nk],np.linspace(.2,.30000001,11))
    old_fsf = saved.LOGMSTAR.to_numpy()-2*np.log10(.704)
    curves = {'GAT':old_curve,'P12-A':new_curve}
    report = {'completeness_curves':curves, 'selection':{}, 'confidence':{}, 'sharpening':{},
              'limits':'Flux-limit mass diagnostic, not proof of total completeness: surface brightness, fibre assignment, SED errors and redshift failures remain. Old FSF completeness tested in its own h=.704 mass units, current in CG15 units; outcome/control masses CG15 throughout. No extrapolation, new model inference or changed tracer selection.',
              'method':'Within each .005 old/.01 current z slice: faintest20% relative to target limit, 90th percentile logM+.4(r-r_limit). Primary comparison main BGS_BRIGHT in identical RA/Dec wedge. Paired cuts use same .1dex CG15mass/.01z cells, reference mixture and hard/soft support; >=20 effective per class.64sky bootstrap draws.'}
    report['selection']['old_main_bright'] = paired_selection(old,old_good,bright & main,old_mass,old_ssfr,'Main BRIGHT versus full old property sample')
    report['selection']['old_common_sky'] = paired_selection(old,old_good & bright & main,footprint,old_mass,old_ssfr,'Same angular wedge versus main BRIGHT all-sky')
    for model,a,base,mass,ssfr,cmass in [('GAT',old,old_base,old_mass,old_ssfr,old_fsf),('P12-A',new,new_base,new_mass,new_ssfr,new_mass)]:
        for extra in [0.,.2]:
            comp = apply_limit(a['z'],cmass,curves[model],extra)
            key = model+('_complete90' if extra==0 else '_complete90_plus0p2')
            report['selection'][key] = paired_selection(a,base,comp,mass,ssfr,key)
        # Confidence selection is a separate intervention on the original sample.
        baseline = old_good if model=='GAT' else new_good
        for cut in [.7,.85,.95]:
            high = np.max(a['P_WEB'],axis=1) >= cut
            report['confidence'][model+'_maxP_'+str(cut)] = paired_selection(a,baseline,high,mass,ssfr,'Confidence cut at '+str(cut))
        ii = np.flatnonzero(baseline); p = a['P_WEB'][ii].astype(float)
        weights = {'p':p,'p_squared':p**2/(p**2).sum(1)[:,None], 'p_fourth':p**4/(p**4).sum(1)[:,None], 'hard':np.eye(4)[p.argmax(1)]}
        cells = np.floor((mass[ii]-9)/.1).astype(int)*100 + np.floor(a['z'][ii]/.01).astype(int)
        region = hp.ang2pix(8,a['RA'][ii],a['DEC'][ii],lonlat=True)
        report['sharpening'][model] = {}
        for name,y in [('ssfr',ssfr[ii]),('quenched',(ssfr[ii]<-11).astype(float))]:
            r,_ = estimate(y,weights,cells,region,np.zeros(len(ii)),nboot=64);r.pop('paired_density_minus_soft_contrast')
            report['sharpening'][model][name] = r
    report['script_sha256'] = sha(__file__)
    dump(OUT/'RESULTS.json',report)
    for category in ['selection','confidence']:
        for name,r in report[category].items():
            values = {k:round(v['contrast']*100,3) for k,v in r.get('properties',{}).get('quenched',{}).get('methods',{}).items()}
            print(category,name,r['status'],values,flush=True)
    print('COMPLETE',OUT,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','analyse']);a=p.parse_args()
    prepare() if a.stage=='prepare' else analyse()
