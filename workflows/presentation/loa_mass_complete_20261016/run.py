"""Parent audit and empirical mass-completeness selection; compute nodes only.

No inference, training, scheduler actions, or overwriting of released catalogues.
Completeness is a flux-limit estimate conditional on observed SED-fit galaxies.
"""
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import socket
import numpy as np

VAC = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')
PHOT = Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
SED = Path('/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits')
OUT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_mass_complete_20261016_v1')
EDGES = np.round(np.arange(.15, .550001, .01), 8)


def write_json(p, value):
    with p.open('x') as f:
        json.dump(value, f, indent=2, allow_nan=False)


def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def limiting_curve(z, mass, r, rlim, sfr, edges=EDGES):
    """95th percentile of rescaled masses in faintest 20%, by population.

Use max(all, quiescent, star-forming) in each bin, with minimum 100 rows
per population. Insufficient support leaves the bin unavailable, not filled.
"""
    valid_sfr = np.isfinite(sfr) & (sfr >= 0)
    q = valid_sfr & (sfr / 10**mass < 1e-11)
    groups = {'all': np.ones(len(z), bool), 'quiescent': q,
              'starforming': valid_sfr & ~q}
    offset = r - rlim
    results = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        base = (z >= lo) & (z < hi) & np.isfinite(mass) & np.isfinite(offset) & (offset < 0)
        item = {'zlo': float(lo), 'zhi': float(hi), 'populations': {}}
        for name, mask in groups.items():
            use = base & mask
            rec = {'rows': int(use.sum()), 'usable': bool(use.sum() >= 100)}
            if rec['usable']:
                faint = use & (offset >= np.quantile(offset[use], .8))
                limits = mass[faint] + .4 * offset[faint]
                rec.update(faint_rows=int(faint.sum()), p90=float(np.quantile(limits, .90)),
                           p95=float(np.quantile(limits, .95)))
            item['populations'][name] = rec
        item['usable'] = all(x['usable'] for x in item['populations'].values())
        item['logmass_limit'] = max(x['p95'] for x in item['populations'].values()) if item['usable'] else None
        results.append(item)
    return results


def selection(z, mass, curve):
    limit = np.full(len(z), np.nan)
    for b in curve:
        if b['usable']:
            limit[(z >= b['zlo']) & (z < b['zhi'])] = b['logmass_limit']
    window = [b for b in curve if b['zlo'] >= .2 and b['zhi'] <= .3]
    if len(window) != 10 or not all(b['usable'] for b in window):
        raise ValueError('Insufficient support for the primary 0.20-0.30 sample')
    fixed = np.ceil(max(b['logmass_limit'] for b in window) * 100) / 100
    return limit, mass >= limit, (z >= .2) & (z < .3) & (mass >= fixed), float(fixed)


def unique_only(a, name='TARGETID'):
    _, inv, count = np.unique(a[name], return_inverse=True, return_counts=True)
    return a[count[inv] == 1], int((count[inv] != 1).sum())


def read_chunks(path, cols, predicate):
    import fitsio
    chunks = []
    with fitsio.FITS(path) as f:
        n = f[1].get_nrows()
        for start in range(0, n, 300000):
            a = f[1].read(rows=np.arange(start, min(start+300000, n)), columns=cols)
            chunks.append(a[predicate(a)])
            if start % 3000000 == 0:
                print(path.name, start, n, flush=True)
    return np.concatenate(chunks)


def prepare(out):
    import fitsio
    from numpy.lib.recfunctions import append_fields
    # Audit the broad SED catalogue first: no environment or mass-threshold cut.
    cols = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'Z', 'SPECTYPE', 'SURVEY',
            'DELTACHI2', 'MASS_CG_15', 'SFR_CG_15']
    cache=out/'sed_parent.npy'
    sed = np.load(cache,allow_pickle=False) if cache.exists() else read_chunks(SED, cols, lambda a: (a['Z'] > 0) & (a['Z'] < .55)
                      & (np.char.strip(a['SPECTYPE'].astype('U')) == 'GALAXY')
                      & (a['DELTACHI2'] > 40)
                      & np.isfinite(a['MASS_CG_15']) & (a['MASS_CG_15'] > 0))
    if not cache.exists():np.save(cache,sed,allow_pickle=False)
    sed, duplicates = unique_only(sed)
    sed.sort(order='TARGETID')
    ze = np.linspace(0, .55, 276); me = np.linspace(6, 13, 281)
    hist = np.histogram2d(sed['Z'], np.log10(sed['MASS_CG_15']), bins=(ze, me))[0]
    low = (sed['Z'] < .15)
    audit = {'sed_unique_valid_mass_rows': len(sed), 'duplicate_rows_excluded': duplicates,
             'low_z_rows': int(low.sum()), 'low_z_below_1e9': int((low & (sed['MASS_CG_15'] < 1e9)).sum()),
             'low_z_1e7_to_1e8': int((low & (sed['MASS_CG_15'] >= 1e7) & (sed['MASS_CG_15'] < 1e8)).sum()),
             'parent_definition': 'Unique IDs among positive-mass GALAXY SED rows at 0<z<.55 and DELTACHI2>40, all targeting programmes; not identical to the historical 5,130,670-object example. SED catalogue lacks ZWARN; Bright parent additionally requires ZWARN=0.'}
    photcols = ['TARGETID', 'RA', 'DEC', 'Z_not4clus', 'FLUX_G', 'FLUX_R',
                'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'PHOTSYS', 'ZWARN', 'DELTACHI2']
    phot = read_chunks(PHOT, photcols, lambda a: (a['Z_not4clus'] > 0) & (a['Z_not4clus'] < .55)
                       & (a['ZWARN']==0) & (a['DELTACHI2']>40))
    phot, pd = unique_only(phot)
    at = np.searchsorted(sed['TARGETID'], phot['TARGETID'])
    found = (at < len(sed)) & (sed['TARGETID'][np.minimum(at, len(sed)-1)] == phot['TARGETID'])
    phot = phot[found]; sedmatch = sed[at[found]]
    # Use the same geometric tolerance as a normal catalogue cross-match;
    # record rather than hide objects that fail the check.
    sep = np.hypot((phot['RA']-sedmatch['TARGET_RA'])*np.cos(np.deg2rad(phot['DEC'])), phot['DEC']-sedmatch['TARGET_DEC'])*3600
    dz = np.abs(phot['Z_not4clus']-sedmatch['Z'])/(1+phot['Z_not4clus'])
    match = (sep <= 1) & (dz <= .001)
    audit.update(bright_duplicate_rows_excluded=pd, bright_id_matches=len(phot),
                 bright_sky_z_rejected=int((~match).sum()))
    phot = phot[match]; sedmatch = sedmatch[match]
    sys = np.char.strip(phot['PHOTSYS'].astype('U')); assert np.isin(sys, ['N','S']).all()
    rlim = np.where(sys == 'N', 19.54, 19.5)
    with np.errstate(invalid='ignore', divide='ignore'):
        r = 22.5-2.5*np.log10(phot['FLUX_R']/phot['MW_TRANSMISSION_R'])
        colour = -2.5*np.log10((phot['FLUX_G']/phot['MW_TRANSMISSION_G'])/(phot['FLUX_R']/phot['MW_TRANSMISSION_R']))
    bright = np.isfinite(r) & (r < rlim) & (phot['MW_TRANSMISSION_R'] > 0)
    phot, sedmatch, r, rlim, colour = [v[bright] for v in (phot, sedmatch, r, rlim, colour)]
    z = phot['Z_not4clus']; mass = np.log10(sedmatch['MASS_CG_15'])
    bh = np.histogram2d(z, mass, bins=(ze, me))[0]
    curve = limiting_curve(z, mass, r, rlim, sedmatch['SFR_CG_15'])
    np.savez(out/'bright_parent.npz',TARGETID=phot['TARGETID'],z=z,mass=mass,r=r,rlim=rlim,sfr=sedmatch['SFR_CG_15'])
    audit.update(bright_valid_mass_rows=len(phot), bright_low_z_below_1e9=int(((z < .15)&(mass < 9)).sum()),
                 bright_low_z_1e7_to_1e8=int(((z < .15)&(mass >= 7)&(mass < 8)).sum()))
    write_json(out/'parent_audit.json', audit)
    write_json(out/'mass_limits.json', curve)
    np.savez_compressed(out/'parent_histograms.npz', all_sed=hist, bright=bh, z_edges=ze, mass_edges=me)
    del sed
    # Read the original VAC without changing its columns or filtering on SFR.
    a = fitsio.read(VAC)
    order = np.argsort(phot['TARGETID']); ids = phot['TARGETID'][order]
    at = np.searchsorted(ids, a['TARGETID']); valid = (at < len(ids)) & (ids[np.minimum(at,len(ids)-1)] == a['TARGETID'])
    ix = order[np.minimum(at,len(ids)-1)]
    valid &= a['SUPPORTED'].astype(bool) & a['CIGALE_MATCHED'].astype(bool) & ~a['CIGALE_AMBIGUOUS'].astype(bool)
    valid &= (np.char.strip(a['CIGALE_SPECTYPE'].astype('U')) == 'GALAXY') & np.isfinite(a['MASS_CG_15']) & (a['MASS_CG_15'] > 0)
    assert np.array_equal(a['Z'][valid], phot['Z_not4clus'][ix[valid]])
    assert np.array_equal(a['RA'][valid], phot['RA'][ix[valid]])
    assert np.array_equal(a['DEC'][valid], phot['DEC'][ix[valid]])
    # A different SED fit for the same ID cannot silently define the selection.
    assert np.array_equal(a['MASS_CG_15'][valid], sedmatch['MASS_CG_15'][ix[valid]])
    rows = np.flatnonzero(valid); a = a[valid]; ix = ix[valid]
    assert len(np.unique(a['TARGETID'])) == len(a)
    mass = np.log10(a['MASS_CG_15']); limits, varying, fixed, threshold = selection(a['Z'], mass, curve)
    sf = a['SFR_CG_15']; sfok = np.isfinite(sf) & (sf >= 0)
    quenched = np.full(len(a), -1, dtype='i2'); quenched[sfok] = (sf[sfok]/a['MASS_CG_15'][sfok] < 1e-11).astype('i2')
    a = append_fields(a, ['SOURCE_ROW','R_DERED','R_LIMIT','GR_OBS','LOGM_LIMIT95','SELECT_ZDEPENDENT','SELECT_FIXED020030','VALID_SFR','QUENCHED_SSFR11'],
                      [rows,r[ix],rlim[ix],colour[ix],limits,varying,fixed,sfok,quenched],usemask=False)
    header = {'MCFIX':threshold, 'MCZLO':.2, 'MCZHI':.3, 'MCQUANT':.95,
              'MCMETH':'max(all,Q,SF) faintest 20%', 'MCCAVEAT':'Empirical flux completeness only'}
    for name, keep in [('zdependent',varying),('fixed020030',fixed)]:
        path = out/f'LOA_P12A_MASS_COMPLETE_{name}.fits'
        if path.exists(): raise FileExistsError(path)
        fitsio.write(path, a[keep], header=header)
    # Small selected input cache, including parent comparison at .20-.30.
    sample = (a['Z'] >= .2) & (a['Z'] < .3)
    np.save(out/'analysis_input.npy', a[sample], allow_pickle=False)
    summary = {'eligible_vac_rows':len(a), 'zdependent_rows':int(varying.sum()),
               'fixed020030_rows':int(fixed.sum()), 'fixed_logmass':threshold,
               'fixed_mass_msun':float(10**threshold), 'invalid_sfr_fixed':int((fixed & ~sfok).sum()),
               'parent020030_rows':int(sample.sum())}
    write_json(out/'selection_summary.json', summary)
    print('SELECTION', summary, flush=True)


def controlled(y, weights, cells, region, minimum=20, nboot=256):
    """Common-stratum standardized means; paired sky-block bootstrap, fixed mix."""
    cell_ids, c = np.unique(cells, return_inverse=True); _, reg = np.unique(region, return_inverse=True)
    nc = c.max()+1; nr = reg.max()+1; flat = c*nr+reg
    C=[]; S=[]; good=np.ones(nc,bool)
    for w in weights.values():
        counts=np.stack([np.bincount(flat,weights=w[:,j],minlength=nc*nr).reshape(nc,nr) for j in range(4)])
        sums=np.stack([np.bincount(flat,weights=w[:,j]*y,minlength=nc*nr).reshape(nc,nr) for j in range(4)])
        w2=np.stack([np.bincount(c,weights=w[:,j]**2,minlength=nc) for j in range(4)])
        good &= (counts.sum(2)**2/np.maximum(w2,1e-30) >= minimum).all(0)
        C.append(counts); S.append(sums)
    if not good.any(): return {'status':'unavailable due to insufficient common mass-redshift support'}
    ref=np.bincount(c,minlength=nc)[good].astype(float);ref/=ref.sum()
    C=np.stack(C)[:,:,good,:]; S=np.stack(S)[:,:,good,:]
    mean=(S.sum(-1)/C.sum(-1)*ref).sum(-1)
    raw_common=S.sum((-1,-2))/C.sum((-1,-2))
    bs=[];rng=np.random.default_rng(20261016)
    for _ in range(nboot):
        mult=np.bincount(rng.integers(nr,size=nr),minlength=nr); den=C@mult
        if (den>0).all():bs.append(((S@mult)/den*ref).sum(-1))
    if len(bs)<.8*nboot:return {'status':'unavailable due to unstable bootstrap support'}
    bs=np.array(bs); d=bs[:,:,-1]-bs[:,:,0]
    return {'status':'measured','rows':int(good[c].sum()),'input_rows':len(y),'cells':int(good.sum()),
            'common_cell_ids':cell_ids[good].tolist(),'reference_weights':ref.tolist(),
            'sky_blocks':int(nr),'bootstrap_valid':len(bs),'methods':{
                name:{'mean':mean[i].tolist(),'raw_on_common_rows':raw_common[i].tolist(),'ci95':np.quantile(bs[:,i],[.025,.975],axis=0).tolist(),
                      'contrast':float(mean[i,-1]-mean[i,0]),'contrast_ci95':np.quantile(d[:,i],[.025,.975]).tolist()}
                for i,name in enumerate(weights)}}


def analyse(out, full=False):
    import healpy as hp
    if full:
        import fitsio
        a=fitsio.read(out/'LOA_P12A_MASS_COMPLETE_zdependent.fits')
    else:a=np.load(out/'analysis_input.npy',allow_pickle=False)
    z=a['Z']; m=np.log10(a['MASS_CG_15']); p=a['P_WEB'].astype(float)
    assert np.isfinite(p).all() and (p>=0).all() and np.allclose(p.sum(1),1,atol=1e-5)
    trace=a['EIGENVALUE_MEAN'].sum(1); assert np.isfinite(trace).all()
    region=hp.ang2pix(8,a['RA'],a['DEC'],lonlat=True)
    sf=a['SFR_CG_15']; sfok=a['VALID_SFR'].astype(bool)
    ss=np.full(len(a),np.nan)
    with np.errstate(divide='ignore',invalid='ignore'): ss[sfok]=np.log10(sf[sfok]/a['MASS_CG_15'][sfok])
    properties={'quenched_ssfr11':(ss < -11).astype(float), 'quenched_ssfr10p5':(ss < -10.5).astype(float),
                'quenched_ssfr11p5':(ss < -11.5).astype(float), 'logmass':m,
                'logssfr_floor14':np.maximum(ss,-14), 'gr_observed':a['GR_OBS'],
                'av':np.where(a['AV_CG_15']>=0,a['AV_CG_15'],np.nan),
                'fagn':np.where((a['FAGN_CG_15']>=0)&(a['FAGN_CG_15']<=1),a['FAGN_CG_15'],np.nan)}
    with np.errstate(divide='ignore',invalid='ignore'):
        properties['logage_Gyr']=np.log10(a['AGE_CG_15']/1000)
        properties['logsfr_floor14']=np.maximum(np.log10(np.where(sfok,sf,np.nan)),-14)
    for name in properties:
        if 'quenched' in name: properties[name][~sfok]=np.nan
    results={}
    selections=[('flux_parent',np.ones(len(a),bool)),('zdependent',a['SELECT_ZDEPENDENT']),('fixed',a['SELECT_FIXED020030']),
                ('fixed_plus0p1',a['SELECT_FIXED020030'] & (m>=json.loads((out/'selection_summary.json').read_text())['fixed_logmass']+.1))]
    if full:selections=[('zdependent_full',np.ones(len(a),bool))]
    for sel, mask in selections:
        results[sel]={}
        for width in [.1,.05]:
            block={}
            for name,y in properties.items():
                use=mask & np.isfinite(y); ii=np.flatnonzero(use)
                zb=np.floor((z[ii]-.2+1e-9)/.01).astype(int)
                cells=zb if name=='logmass' else np.floor((m[ii]-6)/width).astype(int)*100+zb
                # All class curves share the same matched cells. Trace quartiles
                # are an additional descriptor, not equivalent physical classes.
                hard=np.eye(4)[p[ii].argmax(1)]
                weights={'soft_web':p[ii],'hard_web':hard}
                r=controlled(y[ii],weights,cells,region[ii])
                r['raw']={k:(w.T@y[ii]/w.sum(0)).tolist() for k,w in weights.items()}
                r['control']='redshift only' if name=='logmass' else f'{width} dex mass and .01 redshift'
                r['selected_z_range']=[float(z[ii].min()),float(z[ii].max())] if len(ii) else None
                if r['status']=='measured':
                    retained=ii[np.isin(cells,r['common_cell_ids'])]
                    r['matched_z_range']=[float(z[retained].min()),float(z[retained].max())]
                    r['matched_rows_by_zbin']=np.histogram(z[retained],EDGES)[0].tolist()
                block[name]=r
            results[sel][str(width)]=block
            primary=block['quenched_ssfr11']
            print('ANALYSIS',sel,width,primary['status'],primary.get('rows'),
                  {k:(v['contrast'],v['contrast_ci95']) for k,v in primary.get('methods',{}).items()},flush=True)
    write_json(out/('relationships_full.json' if full else 'relationships.json'),results)


def render(out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    d=np.load(out/'parent_histograms.npz'); curves=json.loads((out/'mass_limits.json').read_text())
    summary=json.loads((out/'selection_summary.json').read_text())
    audit=json.loads((out/'parent_audit.json').read_text())
    vmax=max(2,d['all_sed'].max(),d['bright'].max())
    fig,axs=plt.subplots(1,2,figsize=(13,5.5),sharey=True)
    for ax,k,title in zip(axs,['all_sed','bright'],['Broad Loa SED catalogue: before environment cuts','BGS Bright parent: before environment cuts']):
        h=d[k];mesh=ax.pcolormesh(d['z_edges'],d['mass_edges'],np.ma.masked_equal(h.T,0),norm=LogNorm(vmin=1,vmax=vmax),cmap='magma_r',rasterized=True)
        ax.axvline(.15,c='gray',ls='--',label='P12-A lower boundary')
        quality='Positive CG15 mass, Δχ² > 40'+('; ZWARN = 0' if k=='bright' else '')
        ax.set(xlabel='Redshift',title=title+'\n'+quality,xlim=(0,.55),ylim=(6.5,13))
        fig.colorbar(mesh,ax=ax,label='Galaxy count')
        count=audit['sed_unique_valid_mass_rows' if k=='all_sed' else 'bright_valid_mass_rows']
        ax.text(.02,.96,f'N = {count:,}',transform=ax.transAxes,va='top',bbox=dict(fc='white',ec='none',alpha=.8))
    good=[x for x in curves if x['usable']]
    axs[1].plot([(x['zlo']+x['zhi'])/2 for x in good],[x['logmass_limit'] for x in good],c='magenta',label='Conservative empirical 95% limit')
    axs[1].plot([.2,.3],[summary['fixed_logmass']]*2,c='cyan',lw=3,label='Fixed .20-.30 selection')
    axs[0].set_ylabel('log10 stellar mass / solar mass (CG15)');axs[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'parent_mass_redshift.png',dpi=180);fig.savefig(out/'parent_mass_redshift.pdf');plt.close(fig)
    r=json.loads((out/'relationships.json').read_text())
    props=['quenched_ssfr11','logmass','logssfr_floor14','logsfr_floor14','gr_observed','logage_Gyr','av','fagn']
    labels=['Quenched fraction (sSFR < 10⁻¹¹ yr⁻¹)','log₁₀ stellar mass [M☉]',
            'Mean log₁₀ sSFR [yr⁻¹; floor −14]','Mean log₁₀ SFR [M☉/yr; floor −14]',
            'Observed g − r [mag]','log₁₀ fitted age [Gyr]','Fitted dust Aᵥ [mag]','Fitted AGN contribution']
    fig,axs=plt.subplots(2,4,figsize=(17,8))
    for ax,prop,label in zip(axs.flat,props,labels):
        for sel,col,style,legend in [('flux_parent','#666666',':','Flux-limited parent'),('zdependent','#d68a00','--','Redshift-dependent mass cut'),('fixed','#087ca7','-','Fixed mass cut')]:
            rec=r[sel]['0.1'][prop]
            if rec['status']!='measured':continue
            v=rec['methods']['soft_web'];mu=np.array(v['mean']);ci=np.array(v['ci95'])
            ax.plot(range(4),mu,marker='o',ls=style,color=col,label=legend)
            ax.vlines(range(4),ci[0],ci[1],color=col)
        rec=r['fixed']['0.1'][prop]
        note=f"Fixed: {rec['rows']:,} matched" if rec['status']=='measured' else 'Fixed: insufficient support'
        if prop=='logmass':note+=' · z control only'
        ax.set_title(label+'\n'+note,fontsize=9);ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'],rotation=20);ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=8);fig.suptitle('Loa .20 ≤ z < .30: mass/redshift standardized property–environment curves\n95% sky-block intervals; inferred classes; flux-completeness estimate, not causal evidence')
    fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(out/'property_environment.png',dpi=180);fig.savefig(out/'property_environment.pdf');plt.close(fig)
    primary=r['fixed']['0.1']['quenched_ssfr11']
    if primary['status']=='measured':
        fig,axs=plt.subplots(1,2,figsize=(11.5,5))
        for method,col,style in [('soft_web','#087ca7','-'),('hard_web','#d68a00','--')]:
            v=primary['methods'][method];ci=np.array(v['ci95'])*100
            axs[0].plot(range(4),np.array(v['mean'])*100,marker='o',ls=style,c=col,label=method.replace('_',' '))
            axs[0].vlines(range(4),ci[0],ci[1],color=col)
            keys=['quenched_ssfr10p5','quenched_ssfr11','quenched_ssfr11p5']
            vv=[r['fixed']['0.1'][k]['methods'][method] for k in keys]
            off=-.035 if method=='soft_web' else .035;x=np.arange(3)+off
            axs[1].plot(x,[v['contrast']*100 for v in vv],marker='o',ls=style,c=col)
            bounds=np.array([v['contrast_ci95'] for v in vv])*100
            axs[1].vlines(x,bounds[:,0],bounds[:,1],color=col)
        axs[0].set_xticks(range(4),['Void','Wall','Filament','Cluster']);axs[0].set_ylabel('Mass/redshift-standardized quenched fraction [%]')
        axs[0].set_title('Primary: sSFR < 10⁻¹¹ yr⁻¹');axs[0].legend()
        axs[1].set_xticks(range(3),['−10.5','−11 (primary)','−11.5']);axs[1].set_xlabel('log₁₀ sSFR threshold [yr⁻¹]')
        axs[1].set_ylabel('Cluster − void contrast [percentage points]');axs[1].axhline(0,c='gray',lw=1);axs[1].set_title('Quenching-definition sensitivity')
        for ax in axs:ax.grid(alpha=.2)
        fig.suptitle(f"Loa 0.20 ≤ z < 0.30 · log₁₀ M★ ≥ {summary['fixed_logmass']:.2f}\n{primary['rows']:,} common-support galaxies · 95% sky-block intervals")
        fig.tight_layout(rect=[0,0,1,.87]);fig.savefig(out/'quenched_primary.png',dpi=180);fig.savefig(out/'quenched_primary.pdf');plt.close(fig)
    if (out/'relationships_full.json').exists():
        rr=json.loads((out/'relationships_full.json').read_text())['zdependent_full']['0.1']
        fig,axs=plt.subplots(2,4,figsize=(17,8))
        for ax,prop,label in zip(axs.flat,props,labels):
            rec=rr[prop]
            if rec['status']=='measured':
                for method,col,style in [('soft_web','#087ca7','-'),('hard_web','#d68a00','--')]:
                    v=rec['methods'][method];ci=np.array(v['ci95'])
                    ax.plot(range(4),v['mean'],marker='o',ls=style,c=col,label=method.replace('_',' '))
                    ax.vlines(range(4),ci[0],ci[1],color=col)
                note=f"{rec['rows']:,} matched / {rec['input_rows']:,} valid"
                if prop=='logmass':note+=' · z only'
            else:ax.text(.1,.5,rec['status'],transform=ax.transAxes,wrap=True,fontsize=8)
            ax.set_title(label+'\n'+(note if rec['status']=='measured' else 'Insufficient support'),fontsize=9);ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'],rotation=20);ax.grid(alpha=.2)
        axs[0,0].legend(fontsize=8)
        fig.suptitle('Loa selection: 0.15 ≤ z < 0.55 · matched comparison: 0.15 ≤ z < 0.49\nHigher-z cells lack common class support · 95% sky-block intervals')
        fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(out/'property_environment_full.png',dpi=180);fig.savefig(out/'property_environment_full.pdf');plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['prepare','analyse','analyse_full','render','all'],default='all');ap.add_argument('--out',type=Path,default=OUT);args=ap.parse_args()
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Payload requires an allocated compute node')
    args.out.mkdir(parents=True,exist_ok=True)
    if args.stage in ['prepare','all']:prepare(args.out)
    if args.stage in ['analyse','all']:analyse(args.out)
    if args.stage in ['analyse_full','all']:analyse(args.out,full=True)
    if args.stage in ['render','all']:render(args.out)
    if args.stage=='all':
        write_json(args.out/'MANIFEST.json',{'script_sha256':sha(Path(__file__)), 'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),
            'sources':{str(p):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in (VAC,PHOT,SED)},
            'outputs':{p.name:sha(p) for p in args.out.iterdir() if p.is_file()},
            'limitations':'Empirical observed-sample flux completeness, no absolute completeness guarantee; fixed reference sky bootstrap excludes SED and model-transfer systematics. No inference retraining.'})

if __name__=='__main__':main()
