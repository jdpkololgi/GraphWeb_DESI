"""Selection-only z=0 experiment; no lightcone, evolution, or new inference."""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ[key] = '1'
from pathlib import Path
import json
import hashlib
import csv
import numpy as np
import h5py
from astropy.cosmology import Planck18
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from loa_tng_abacus_overlay_20261015 import property_data

ROOT = Path('/global/u2/d/dkololgi/fmsu_tokyo_20261010')
T = ROOT / 'tng_density_control_20261015'
OUT = ROOT / 'tng_bright_selection_20261015'
OUT.mkdir(exist_ok=True)
CAT = Path('/pscratch/sd/d/dkololgi/TNG300-1/groups_099')

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def quarters(x):
    _, inv, n = np.unique(x, return_inverse=True, return_counts=True)
    end = np.cumsum(n); start = end - n; size = len(x)
    return np.array([np.maximum(0, np.minimum(end, (k+1)*size/4)
                                - np.maximum(start, k*size/4))/n
                     for k in range(4)]).T[inv]

def main():
    d = np.load(T/'catalogue.npz'); f = np.load(T/'features.npz')
    ids = f['test_ids']; mass = np.log10(d['mass_msun'][ids])
    ids = ids[(mass >= 9) & (mass < 11.5) & d['valid_query'][ids]]
    sid = d['subhalo_id'][ids]
    cache = OUT/'photometry.npz'
    if not cache.exists():
        # Read only the photometric columns needed, joining by global Subfind ID.
        phot = np.full((len(ids), 2), np.nan); offset = 0; sources = []
        for part in range(600):
            p = CAT/f'fof_subhalo_tab_099.{part}.hdf5'
            with h5py.File(p, 'r') as h:
                n = int(h['Header'].attrs['Nsubgroups_ThisFile'])
                assert abs(float(h['Header'].attrs['Redshift'])) < 1e-8
                ii = np.flatnonzero((sid >= offset) & (sid < offset+n))
                if len(ii):
                    values = h['Subhalo/SubhaloStellarPhotometrics'][:, 4:6]
                    phot[ii] = values[sid[ii]-offset]
                offset += n
                sources.append({'path':str(p), 'bytes':p.stat().st_size,
                                'mtime_ns':p.stat().st_mtime_ns, 'rows':len(ii)})
            if part % 100 == 0: print('photometry chunk', part, flush=True)
        assert np.isfinite(phot).all()
        np.savez_compressed(cache, ids=ids, subhalo_id=sid, g=phot[:,0], r=phot[:,1])
        (OUT/'photometry_sources.json').write_text(json.dumps(sources, indent=2))
    phot = np.load(cache)
    assert np.array_equal(phot['ids'], ids) and np.array_equal(phot['subhalo_id'], sid)
    assert np.allclose(phot['g']-phot['r'], d['g_r'][ids], atol=2e-6)
    mass = np.log10(d['mass_msun'][ids]); mb = np.floor((mass-9)/.25).astype(int)
    q = (d['sfr'][ids]/d['mass_msun'][ids] < 1e-11).astype(float)
    trace = np.load(T/'truth.npz')['eig_R7_n512'][ids].sum(1)
    block = np.ravel_multi_index(np.minimum((d['pos_hmpc'][ids]/205*5).astype(int),4).T,(5,5,5))
    a = np.load(ROOT/'loa_tng_abacus_overlay_20261015/loa_predictions.npz')
    lp = np.load(ROOT/'legacy_selection_tests_20261014/p12a_targeting.npz')
    assert np.array_equal(a['TARGETID'], lp['TARGETID'])
    good,y = property_data(a)
    good &= (lp['r'] < lp['rlimit']) & (y['mass'] < 11.5) & (a['z'] >= .2) & (a['z'] < .3)
    assert np.all(lp['rlimit'][good] == 19.5), 'This implementation expects southern wedge.'
    lm = np.floor((y['mass'][good]-9)/.25+1e-10).astype(int)
    lz = np.floor((a['z'][good]-.2)/.01+1e-8).astype(int)
    lq = y['quenched'][good]; lt = a['abacus'][good].sum(1)
    _, lb = np.unique(a['region'][good], return_inverse=True)
    nb = lb.max()+1
    td = np.random.default_rng(731).multinomial(125,np.ones(125)/125,size=200)
    ld = np.random.default_rng(732).multinomial(nb,np.ones(nb)/nb,size=200)
    checks = 0
    def stat(yy, xx, bb, draws, weights=None):
        nonlocal checks
        w = quarters(xx) if weights is None else weights
        nblocks = draws.shape[1]
        den = np.array([np.bincount(bb, weights=w[:,k], minlength=nblocks) for k in range(4)]).T
        num = np.array([np.bincount(bb, weights=w[:,k]*yy, minlength=nblocks) for k in range(4)]).T
        val = num.sum(0)/den.sum(0)
        assert np.allclose(val, (w.T@yy)/w.sum(0), atol=1e-12)
        db = draws@den
        assert np.all(db>0), 'Unsupported bootstrap cell'
        checks += 1
        return {'n':len(yy), 'mean':val, 'boot':(draws@num)/db}
    arms = ['uncut','bright_color','bright_no_k','bright_color_plus03']
    stats = {k:{} for k in arms+['loa']}; fixed = {}; cells=[]
    for zbin in range(10):
        z = .205 + .01*zbin
        dm = float(Planck18.distmod(z).value)
        # Approximate f_nu power law between SDSS g and r, evaluated at
        # lambda_r/(1+z). Broad filters and spectral breaks are not integrated.
        kcorr = (phot['g']-phot['r'])*np.log(1+z)/np.log(6231/4770)-2.5*np.log10(1+z)
        r0 = phot['r']+dm
        masks = {'uncut':np.ones(len(ids),bool), 'bright_color':r0+kcorr<19.5,
                 'bright_no_k':r0<19.5, 'bright_color_plus03':r0+kcorr+.3<19.5}
        assert np.all(~masks['bright_color_plus03'] | masks['bright_color'])
        for b in range(10):
            li = np.flatnonzero((lm==b)&(lz==zbin)); parent = np.flatnonzero(mb==b)
            indices = {arm:np.flatnonzero((mb==b)&mask) for arm,mask in masks.items()}
            counts = {arm:len(ii) for arm,ii in indices.items()}
            key=(b,zbin); supported = len(li)>=80 and min(counts.values())>=80
            cells.append({'mass_bin':b,'z_bin':zbin,'z_mid':z,'loa':len(li),**counts,'supported':supported})
            if not supported: continue
            stats['loa'][key] = stat(lq[li], lt[li], lb[li], ld)
            for arm,ii in indices.items():
                stats[arm][key] = stat(q[ii],trace[ii],block[ii],td)
            parent_w = quarters(trace[parent]); keep = masks['bright_color'][parent]
            # Same parent rank boundaries: diagnostic of rank redefinition.
            if np.all(parent_w[keep].sum(0)>20):
                fixed[key] = stat(q[parent][keep],trace[parent][keep],block[parent][keep],td,parent_w[keep])
    results={}; table=[]
    def aggregate(s, keys, weights):
        v=sum(w*s[k]['mean'] for k,w in zip(keys,weights))
        boot=sum(w*s[k]['boot'] for k,w in zip(keys,weights))
        return v,boot
    fig,axes=plt.subplots(1,2,figsize=(12,5.9),sharey=True)
    colors=['#333333','#007f98','#888888','#c57613','#a42289']
    labels=['TNG: mass–redshift matched, no brightness cut','TNG: approximate Bright cut',
            'TNG: Bright cut, no K correction','TNG: Bright cut + 0.3 mag dimming','Loa: P12-A ranks, CG15 SFR']
    styles=['o--','s-','.:','^:','D-']
    for panel,high in enumerate([False,True]):
        keys=sorted(k for k in stats['loa'] if not high or k[0]>=6)
        assert keys
        weights=np.array([stats['loa'][k]['n'] for k in keys],float); weights/=weights.sum()
        curves={}; boots={}
        for arm in arms+['loa']:
            v,bt=aggregate(stats[arm],keys,weights);boots[arm]=bt
            curves[arm]={'mean':v.tolist(),'ci95':np.quantile(bt,[.025,.975],axis=0).tolist(),
                         'contrast_pp':float(100*(v[-1]-v[0]))}
            for quarter in range(4):table.append({'panel':'high' if high else 'full','arm':arm,'quarter':quarter+1,'fraction':v[quarter]})
        for arm in arms:
            gap=np.array(curves[arm]['mean'])-curves['loa']['mean']
            curves[arm]['gap_pp']=(100*gap).tolist()
            curves[arm]['mean_gap_pp']=float(100*gap.mean())
            curves[arm]['mean_gap_ci95_pp']=(100*np.quantile((boots[arm]-boots['loa']).mean(1),[.025,.975])).tolist()
        reduction=100*(boots['uncut']-boots['bright_color']).mean(1)
        fixed_keys=[k for k in keys if k in fixed]
        fixed_result=None
        if len(fixed_keys)==len(keys):
            v,bt=aggregate(fixed,keys,weights)
            fixed_result={'mean':v.tolist(),'ci95':np.quantile(bt,[.025,.975],axis=0).tolist()}
        results['high' if high else 'full']={'curves':curves,'cells':[list(k) for k in keys],
            'cell_weights':weights.tolist(),'loa_rows':sum(stats['loa'][k]['n'] for k in keys),
            'loa_eligible':int(np.sum(lm>=6)) if high else len(lm),
            'mean_gap_reduction_pp':float(np.mean(curves['uncut']['gap_pp'])-np.mean(curves['bright_color']['gap_pp'])),
            'paired_reduction_ci95_pp':np.quantile(reduction,[.025,.975]).tolist(),
            'fixed_parent_quarters_bright':fixed_result}
        ax=axes[panel]
        for arm,label,color,style in zip(arms+['loa'],labels,colors,styles):
            c=curves[arm];v=np.array(c['mean'])*100;ci=np.array(c['ci95'])*100
            ax.plot(range(1,5),v,style,color=color,label=label,lw=1.8)
            if arm in ['uncut','bright_color','loa']:ax.fill_between(range(1,5),ci[0],ci[1],color=color,alpha=.10)
        ax.set_title(('All supported masses' if not high else 'High mass: 10.5 ≤ log M < 11.5')+f'\n{len(keys)} shared cells; {results["high" if high else "full"]["loa_rows"]:,} Loa galaxies',fontsize=11)
        ax.set_xticks([1,2,3,4],['Lowest','2','3','Highest']);ax.grid(alpha=.15)
        ax.set_xlabel('Environment quarter within mass–redshift cells')
    axes[0].set_ylabel('Quenched fraction (%)');axes[0].set_ylim(0,100)
    fig.suptitle('Does approximate BGS brightness selection close the TNG–Loa offset?',fontsize=14)
    fig.legend(*axes[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.93),ncol=2,fontsize=9,frameon=False)
    fig.text(.5,.02,'Selection-only experiment: TNG z = 0 placed at observing z = 0.20–0.30; r < 19.5; dust-free g/r photometry.\nAll curves use identical Loa weights in 0.25-dex mass × 0.01-redshift cells (≥80 rows in every arm).\nTNG R7 truth versus Loa P12-A posterior-mean trace ranks; instantaneous versus CG15 SFR remains unmatched.\nShading: conditional 95% spatial bootstrap; repeated TNG placements share the same cube resampling.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.17,1,.77]);fig.savefig(OUT/'selection_overlay.png',dpi=180);fig.savefig(OUT/'selection_overlay.pdf');plt.close(fig)
    receipt={'results':results,'cells':cells,'direct_mean_checks':checks,'tng_unique_parent':len(ids),
        'photometry_color_parity':True,'selection_nested_check':True,'bootstrap_draws':200,
        'method':'Same common support across uncut, approximate color K, zero K, +0.3mag and Loa. Cell midpoints for observing z; same Loa cell weights for every arm; quarters reranked in each arm; fixed-parent diagnostic also saved. No independent copies, no evolution, dust, noise, surface-brightness/fiber selection or graph rebuilding. Dimming is sensitivity only, not a dust model.',
        'photometry_reference':'https://www.tng-project.org/data/docs/specifications/',
        'sources':{str(p):sha(p) for p in [T/'catalogue.npz',T/'truth.npz',T/'features.npz',cache,ROOT/'loa_tng_abacus_overlay_20261015/loa_predictions.npz',ROOT/'legacy_selection_tests_20261014/p12a_targeting.npz',Path(__file__)]}}
    (OUT/'RESULTS.json').write_text(json.dumps(receipt,indent=2,allow_nan=False))
    with (OUT/'curves.csv').open('w') as h:
        writer=csv.DictWriter(h,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    print(json.dumps(results,indent=2),flush=True)

if __name__=='__main__':main()
