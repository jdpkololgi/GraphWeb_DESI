"""Joint mass/redshift standardised observed-colour histograms, fixed sky bootstrap."""
import os
os.environ['MPLBACKEND']='Agg'
import json
import hashlib
from pathlib import Path
import socket
from concurrent.futures import ProcessPoolExecutor
import numpy as np

ROOT=Path(__file__).resolve().parents[3]
SOURCE=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_rest_colour_20261016_v1')
OUT=SOURCE.parent/'loa_matched_colour_20261016_v1'
WINDOWS=[('z015020',.15,.20,10.60,534070),('z020022',.20,.22,10.69,232024),('z020030',.20,.30,10.96,611197)]
EDGES=np.linspace(-.5,2.5,91)
NB=len(EDGES)+1  # Explicit underflow/overflow: all finite colours enter normalization.
ENV=['Void','Wall','Filament','Cluster']
COLORS=['#0072B2','#E69F00','#009E73','#CC79A7']
STYLES=['-','--','-.',':']

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()

def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def common_support(m,z,p,dm,dz,minimum=20):
    # Absolute grid edges keep the same cell definition between windows.
    cells=np.column_stack((np.floor(m/dm+1e-9).astype(int),np.floor(z/dz+1e-9).astype(int)))
    ids,c=np.unique(cells,axis=0,return_inverse=True)
    hard=np.eye(4)[p.argmax(1)];good=np.ones(len(ids),bool)
    for w in [hard,p]:
        sums=np.array([np.bincount(c,weights=w[:,j],minlength=len(ids)) for j in range(4)])
        squares=np.array([np.bincount(c,weights=w[:,j]**2,minlength=len(ids)) for j in range(4)])
        good &= (sums**2/np.maximum(squares,1e-30)>=minimum).all(0)
    if not good.any():raise ValueError('No common support')
    keep=good[c]; _, cc=np.unique(c[keep],return_inverse=True)
    target=np.bincount(cc).astype(float);target/=target.sum()
    return keep,cc,ids[good],target

def estimate(y,m,z,w,c,region,target,multipliers):
    nc=len(target);nr=multipliers.shape[1]
    b=np.searchsorted(EDGES,y,side='right')
    raw=np.array([np.bincount(b,weights=w[:,j],minlength=NB)/w[:,j].sum() for j in range(4)])
    probs=[];boots=[];means=[];bootmeans=[];balance=[]
    for j in range(4):
        ww=w[:,j];counts=np.bincount(c,weights=ww,minlength=nc)
        weights=ww*target[c]/counts[c]
        achieved=np.bincount(c,weights=weights,minlength=nc)
        assert np.allclose(achieved,target,atol=1e-12,rtol=1e-10)
        hist=np.bincount(b,weights=weights,minlength=NB)
        assert np.isclose(hist.sum(),1,atol=1e-12)
        # region x cell x colour, including tails; independently reproduces weighted histogram.
        tensor=np.bincount((region*nc+c)*NB+b,weights=ww,minlength=nr*nc*NB).reshape(nr,nc,NB)
        den=tensor.sum((0,2))
        independent=(tensor.sum(0)/den[:,None]*target[:,None]).sum(0)
        assert np.allclose(hist,independent,atol=1e-12)
        boot=(multipliers@tensor.reshape(nr,-1)).reshape(len(multipliers),nc,NB)
        bootden=boot.sum(2);valid=(bootden>0).all(1)
        bootprob=np.full((len(multipliers),NB),np.nan)
        bootprob[valid]=(boot[valid]/bootden[valid,:,None]*target[None,:,None]).sum(1)
        moment=np.bincount(region*nc+c,weights=ww*y,minlength=nr*nc).reshape(nr,nc)
        bm=np.full(len(multipliers),np.nan)
        bm[valid]=((multipliers[valid]@moment)/bootden[valid]*target).sum(1)
        probs.append(hist);boots.append(bootprob);means.append(float(weights@y));bootmeans.append(bm)
        balance.append(dict(neff=float(1/(weights@weights)),max_weight=float(weights.max()),
            mean_logmass=float(weights@m),mean_z=float(weights@z),
            max_cell_error=float(np.max(np.abs(achieved-target))),achieved_cell_weights=achieved.tolist()))
    probs=np.array(probs);boots=np.stack(boots,axis=1);bootmeans=np.stack(bootmeans,axis=1)
    valid=np.isfinite(boots).all((1,2)) & np.isfinite(bootmeans).all(1)
    assert valid.sum()>=205, 'Fewer than 80% valid sky bootstrap draws'
    boots=boots[valid];bootmeans=bootmeans[valid]
    assert np.allclose(boots.sum(2),1,atol=1e-12)
    residual=probs-probs.mean(0)
    bres=boots-boots.mean(1,keepdims=True)
    delta=bootmeans[:,3]-bootmeans[:,0]
    return dict(raw_common=raw.tolist(),matched=probs.tolist(),residual=residual.tolist(),
        residual_ci95=np.quantile(bres,[.025,.975],axis=0).tolist(),
        mean_colour=means,cluster_minus_void=float(means[3]-means[0]),
        cluster_minus_void_ci95=np.quantile(delta,[.025,.975]).tolist(),
        bootstrap_valid=int(valid.sum()),balance=balance,
        max_abs_cdf_vs_pool=np.max(np.abs(np.cumsum(residual,axis=1)),axis=1).tolist())

def analyse(window):
    import fitsio
    import healpy as hp
    key,lo,hi,lim,expected=window
    path=SOURCE/f'observed_original_{key}.fits'
    a=fitsio.read(path,columns=['TARGETID','MASS_CG_15','Z','GR_OBS','P_WEB','RA','DEC'])
    assert len(a)==expected and len(np.unique(a['TARGETID']))==expected
    m=np.log10(a['MASS_CG_15']);z=a['Z'];y=a['GR_OBS'];p=a['P_WEB'].astype(float)
    assert np.isfinite(m).all() and np.isfinite(z).all() and np.isfinite(y).all()
    assert (m>=lim-1e-8).all() and (z>=lo).all() and (z<hi).all()
    assert np.isfinite(p).all() and (p>=0).all() and np.allclose(p.sum(1),1,atol=1e-5)
    hard=np.eye(4)[p.argmax(1)]
    result=dict(window=key,zlo=lo,zhi=hi,mass_limit=lim,rows=len(a),source=str(path),source_sha256=sha(path),grids={})
    for grid,dm,dz in [('primary',.1,.01),('fine',.05,.005)]:
        keep,c,ids,target=common_support(m,z,p,dm,dz)
        blocks,region=np.unique(hp.ang2pix(8,a['RA'][keep],a['DEC'][keep],lonlat=True),return_inverse=True)
        rng=np.random.default_rng(20261016)
        mult=np.array([np.bincount(rng.integers(len(blocks),size=len(blocks)),minlength=len(blocks)) for _ in range(256)],dtype=float)
        rec=dict(dm=dm,dz=dz,rows=int(keep.sum()),retained_fraction=float(keep.mean()),cells=len(ids),
            cell_ids=ids.tolist(),reference_weights=target.tolist(),sky_blocks=len(blocks),
            hard_rows=np.bincount(p[keep].argmax(1),minlength=4).tolist(),methods={})
        for method,w in [('hard',hard),('soft',p)]:
            ans=estimate(y[keep],m[keep],z[keep],w[keep],c,region,target,mult)
            b=np.searchsorted(EDGES,y,side='right')
            ans['raw_all']=[(np.bincount(b,weights=w[:,j],minlength=NB)/w[:,j].sum()).tolist() for j in range(4)]
            rec['methods'][method]=ans
            print(key,grid,method,'rows',rec['rows'],'delta',ans['cluster_minus_void'],'CI',ans['cluster_minus_void_ci95'],flush=True)
        result['grids'][grid]=rec
    write(OUT/f'{key}.json',result)
    return key

def render(figdir_name="figures", csv_name="histograms.csv"):
    import csv
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    figdir=OUT/figdir_name;figdir.mkdir(exist_ok=True)
    records=[json.loads((OUT/f'{w[0]}.json').read_text()) for w in WINDOWS]
    widths=np.diff(EDGES);centres=(EDGES[:-1]+EDGES[1:])/2
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    rows=[]
    for grid in ['primary','fine']:
        for method in ['hard','soft']:
            fig,axs=plt.subplots(3,3,figsize=(15,11),sharex=True,sharey='row')
            for k,r in enumerate(records):
                g=r['grids'][grid];d=g['methods'][method]
                raw=np.array(d['raw_common']);matched=np.array(d['matched']);res=np.array(d['residual']);ci=np.array(d['residual_ci95'])
                for j in range(4):
                    for row,arr in [(0,raw),(1,matched)]:
                        axs[row,k].stairs(arr[j,1:-1]/widths,EDGES,color=COLORS[j],ls=STYLES[j],lw=1.7)
                    axs[2,k].stairs(res[j,1:-1]/widths,EDGES,color=COLORS[j],ls=STYLES[j],lw=1.4)
                    axs[2,k].fill_between(centres,ci[0,j,1:-1]/widths,ci[1,j,1:-1]/widths,color=COLORS[j],alpha=.12,step='mid')
                    for b in range(len(widths)):
                        rows.append(dict(window=r['window'],grid=grid,method=method,environment=ENV[j],
                            colour_lo=float(EDGES[b]),colour_hi=float(EDGES[b+1]),raw_common_density=raw[j,b+1]/widths[b],
                            matched_density=matched[j,b+1]/widths[b],residual_density=res[j,b+1]/widths[b],
                            residual_ci_low=ci[0,j,b+1]/widths[b],residual_ci_high=ci[1,j,b+1]/widths[b]))
                axs[0,k].set_title(f"{r['zlo']:.2f} ≤ z < {r['zhi']:.2f} · log M★ ≥ {r['mass_limit']:.2f}\nCommon support: {g['rows']:,} ({g['retained_fraction']:.1%}); {g['cells']} cells",fontsize=11)
                delta=d['cluster_minus_void'];low,high=d['cluster_minus_void_ci95']
                axs[1,k].text(.03,.96,f'Cluster − void mean colour\n{delta:+.4f} [{low:+.4f}, {high:+.4f}] mag',transform=axs[1,k].transAxes,va='top',fontsize=9)
                axs[2,k].axhline(0,color='#444444',lw=.8)
                axs[2,k].set_xlabel('Observed g − r [mag]')
                for ax in axs[:,k]:
                    ax.set_xlim(.3,1.9);ax.grid(axis='y',alpha=.18)
                axs[0,k].set_ylim(bottom=0);axs[1,k].set_ylim(bottom=0)
            for row, field in [(0, 'raw_common'), (1, 'matched')]:
                maximum=max(np.max(np.array(r['grids'][grid]['methods'][method][field])[:,1:-1]/widths) for r in records)
                axs[row,0].set_ylim(0, maximum*1.26)
            residual_max=max(np.max(np.abs(np.array(r['grids'][grid]['methods'][method]['residual_ci95'])[:,:,1:-1]/widths)) for r in records)
            axs[2,0].set_ylim(-residual_max*1.08,residual_max*1.08)
            axs[0,0].set_ylabel('Before weighting\nCommon-support density [mag⁻¹]')
            axs[1,0].set_ylabel('Mass/redshift matched\nDensity [mag⁻¹]')
            axs[2,0].set_ylabel('Matched − pooled matched\nDensity difference [mag⁻¹]')
            handles=[Line2D([],[],color=c,ls=s,lw=2,label=e) for c,s,e in zip(COLORS,STYLES,ENV)]
            fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.925),ncol=4,frameon=False)
            mode='Highest-probability classes' if method=='hard' else 'Probability-weighted classes'
            dm=records[0]['grids'][grid]['dm'];dz=records[0]['grids'][grid]['dz']
            fig.suptitle(f'Observed colour after joint stellar-mass/redshift matching\n{mode} · cells {dm:.2f} dex × {dz:.3f} in z',fontsize=15,y=.99)
            fig.subplots_adjust(left=.08,right=.99,top=.83,bottom=.115,wspace=.15,hspace=.17)
            fig.text(.08,.025,'Same retained galaxies before/after weighting; all four environments share the pooled joint cell distribution.\n'
                'Bands: pointwise conditional 95% sky-block bootstrap, refitting weights with fixed cells/target. Pooled curve: equal mean of four classes.\n'
                'All finite colours enter normalization; display zoom 0.3–1.9 mag. Observed colours/CIGALE masses; no causal or calibration claim.',fontsize=10)
            for ext in ['png','pdf']:fig.savefig(figdir/f'matched_colour_{grid}_{method}.{ext}',dpi=170)
            plt.close(fig)
    with (OUT/csv_name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Approved compute node required')
    OUT.mkdir(exist_ok=False)
    with ProcessPoolExecutor(max_workers=3) as pool:list(pool.map(analyse,WINDOWS))
    render()
    write(OUT/'MANIFEST.json',dict(job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),script_sha256=sha(Path(__file__)),
        files={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file()}))
    print('COMPLETE',OUT,flush=True)

if __name__=='__main__':main()
