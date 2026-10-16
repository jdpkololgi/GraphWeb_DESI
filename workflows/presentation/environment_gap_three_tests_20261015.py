"""Bounded, single-thread complementary diagnostics; no inference or scheduler changes."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
from pathlib import Path
import json, hashlib
import numpy as np
from scipy.stats import rankdata
import fitsio, h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HOME=Path('/global/u2/d/dkololgi'); ROOT=HOME/'fmsu_tokyo_20261010'
OUT=HOME/'GraphWeb_DESI/docs/figures/environment_gap_three_tests_20261015'; OUT.mkdir(exist_ok=True)
def dump(n,x): (OUT/n).write_text(json.dumps(x,indent=2,allow_nan=False))
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
RNG=np.random.default_rng(20261015)

def test1():
    root=ROOT/'tng_density_control_20261015'
    d=np.load(root/'catalogue.npz'); t=np.load(root/'truth.npz'); f=np.load(root/'features.npz')
    ids=f['test_ids']; mass=np.log10(d['mass_msun'][ids]); keep=(mass>=9)&(mass<11.5)
    ids=ids[keep]; mass=mass[keep]; bins=np.floor((mass-9)/.25).astype(int)
    y=np.column_stack([d['sfr'][ids]/d['mass_msun'][ids]<1e-11,d['g_r'][ids],np.log10(np.maximum(d['sfr'][ids]/d['mass_msun'][ids],1e-14))])
    block=np.ravel_multi_index(np.minimum((d['pos_hmpc'][ids]/205*5).astype(int),4).T,(5,5,5)); _,block=np.unique(block,return_inverse=True); nb=block.max()+1
    bw=RNG.multinomial(nb,np.ones(nb)/nb,size=400)
    # Central is the first subhalo in each FoF group. This census is only TNG z=0.
    first=[]
    for p in sorted(Path('/pscratch/sd/d/dkololgi/TNG300-1/groups_099').glob('fof_subhalo_tab_099.*.hdf5')):
        with h5py.File(p) as h:
            if 'Group/GroupFirstSub' in h: first.append(h['Group/GroupFirstSub'][:])
    central=np.isin(d['subhalo_id'][ids],np.concatenate(first)); satellite=~central
    result={'rows':len(ids),'properties':['quenched_fraction','g_r','log_ssfr'],'arms':[], 'source_hashes':{n:sha(root/n) for n in ['catalogue.npz','truth.npz','features.npz']}}
    fig,axes=plt.subplots(1,3,figsize=(13,4))
    for radius,color in [(2,'#167d9a'),(7,'#d66e27')]:
        eig=t[f'eig_R{radius}_n512'][ids]; trace=eig.sum(1); cls=(eig>.2).sum(1)
        membership=np.zeros(len(ids),int)
        for b in np.unique(bins):
            sel=bins==b; q=np.quantile(trace[sel],[.25,.75]); membership[sel&(trace<=q[0])]=-1;membership[sel&(trace>=q[1])]=1
        for high in [False,True]:
            use=mass>=10.5 if high else np.ones(len(ids),bool)
            for estimator in ['mass_controlled_trace_quartiles','raw_trace_quartiles','mass_controlled_physical_classes']:
                cells=np.unique(bins[use]) if estimator!='raw_trace_quartiles' else [None]
                acc=[]; boots=[]; counts=[]
                for b in cells:
                    cell=use if b is None else use&(bins==b)
                    lo=cell&(cls==0) if 'physical' in estimator else cell&(membership==-1)
                    hi=cell&(cls==3) if 'physical' in estimator else cell&(membership==1)
                    if min(lo.sum(),hi.sum())<20:continue
                    vals=[]; bs=[]
                    for sel in [lo,hi]:
                        c=np.bincount(block[sel],minlength=nb); sums=np.array([np.bincount(block[sel],weights=y[sel,j],minlength=nb) for j in range(3)]).T
                        vals.append(sums.sum(0)/c.sum()); den=bw@c; bs.append(np.divide(bw@sums,den[:,None],out=np.full((400,3),np.nan),where=den[:,None]>0))
                    acc.append(vals[1]-vals[0]);boots.append(bs[1]-bs[0]); counts.append(cell.sum())
                weights=np.array(counts)/sum(counts); point=weights@np.array(acc);boot=np.einsum('i,ijk->jk',weights,np.array(boots))
                record=dict(radius=radius,high_mass=high,estimator=estimator,n=int(use.sum()),support_rows=int(sum(counts)),contrast=point.tolist(),ci95=np.nanquantile(boot,[.025,.975],axis=0).T.tolist(),satellite_fraction=float(satellite[use].mean()),quenched_fraction=float(y[use,0].mean()))
                result['arms'].append(record)
                if estimator=='mass_controlled_trace_quartiles':
                    for j,ax in enumerate(axes):
                        scale=100 if j==0 else 1; x=int(high)+(radius==7)*.17-.085; ci=np.array(record['ci95'][j])*scale
                        ax.errorbar(x,point[j]*scale,yerr=[[point[j]*scale-ci[0]],[ci[1]-point[j]*scale]],fmt='o',capsize=4,color=color,label=f'R{radius}' if not high else None)
    for ax,label in zip(axes,['Quenched contrast (pp)','g − r contrast (mag)','log sSFR contrast (dex)']):
        ax.set(xticks=[0,1],xticklabels=['9 ≤ log M < 11.5','10.5 ≤ log M < 11.5'],ylabel=label);ax.axhline(0,color='gray',lw=.5);ax.legend()
    fig.suptitle('Test 1: true TNG environment, no ML — mass-controlled Q4 minus Q1');fig.tight_layout();fig.savefig(OUT/'test1_ground_truth.png',dpi=160);plt.close(fig)
    dump('test1.json',result);print('Test 1 complete',flush=True)

def extract():
    cache=ROOT/'environment_gap_three_tests_20261015'; cache.mkdir(exist_ok=True)
    target=cache/'uncertainty_sample.npz'
    if target.exists():return np.load(target)
    src=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')
    cols=['TARGETID','RA','DEC','Z','SUPPORTED','CIGALE_MATCHED','CIGALE_AMBIGUOUS','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15','P_WEB','EIGENVALUE_MEAN','EIGENVALUE_Q16','EIGENVALUE_Q84','BOUNDARY_MPC','QUALITY']
    chunks=[];rows=[]
    with fitsio.FITS(src) as f:
        n=f[1].get_nrows()
        for start in range(0,n,100000):
            a=f[1].read(rows=np.arange(start,min(start+100000,n)),columns=cols)
            keep=(a['RA']>=120)&(a['RA']<160)&(a['DEC']>=14.5)&(a['DEC']<30.6)&(a['Z']>=.15)&(a['Z']<.50)&a['SUPPORTED'].astype(bool)
            chunks.append(a[keep]); rows.append(np.flatnonzero(keep)+start)
    a=np.concatenate(chunks);np.savez_compressed(target,**{k:a[k] for k in cols},source_row=np.concatenate(rows))
    dump('extraction.json',dict(source=str(src),source_rows=n,extracted_rows=len(a),source_size=src.stat().st_size,source_mtime_ns=src.stat().st_mtime_ns,selection='SUPPORTED; RA[120,160), DEC[14.5,30.6), z[.15,.50)',extracted_sha256=sha(target),full_source_hash_recomputed=False))
    print('Extracted',len(a),'rows',flush=True);return np.load(target)

def test3():
    d=extract(); m=d['MASS_CG_15']; sfr=d['SFR_CG_15']; p=d['P_WEB']; eig=d['EIGENVALUE_MEAN']; width=(d['EIGENVALUE_Q84'][:,0]-d['EIGENVALUE_Q16'][:,0])/2
    # CG15 MASS is linear solar masses; SFR is linear Msun/year.
    mass=np.log10(np.maximum(m,1e-30)); z=d['Z']; q=(sfr/np.maximum(m,1e-30)<1e-11).astype(float)
    valid=d['CIGALE_MATCHED'].astype(bool)&~d['CIGALE_AMBIGUOUS'].astype(bool)&(np.char.strip(d['CIGALE_SPECTYPE'].astype(str))=='GALAXY')&(m>0)&(sfr>=0)&np.isfinite(m)&np.isfinite(sfr)&(mass>=9)&(mass<12)&np.isfinite(eig).all(1)&np.isfinite(width)&(width>0)&np.isfinite(p).all(1)&(p>=0).all(1)&(np.abs(p.sum(1)-1)<1e-4)
    cells=np.floor((np.where(np.isfinite(mass),mass,0)-9)/.2+1e-7).astype(int)+100*np.floor((z-.15)/.05+1e-7).astype(int)
    low=np.zeros(len(z),bool)
    for c in np.unique(cells[valid]):
        ii=np.flatnonzero(valid&(cells==c)); order=np.lexsort((d['TARGETID'][ii],width[ii]));low[ii[order[:len(ii)//2]]]=True
    assert len(np.unique(d['TARGETID']))==len(z)
    blocks=(np.floor((d['RA']-120)/5).astype(int)+8*np.floor((d['DEC']-14.5)/4).astype(int));_,blocks=np.unique(blocks,return_inverse=True);nb=blocks.max()+1;bw=RNG.multinomial(nb,np.ones(nb)/nb,size=200)
    result=dict(valid_rows=int(valid.sum()),low_width_rows=int((valid&low).sum()),metric='half central 68% interval width of lambda1; NOT posterior standard deviation',selection='lowest exactly floor(Ncell/2) width within 0.2 dex/0.05 z cell; ties by TARGETID',windows=[])
    windows=[(.15,.50),(.15,.25)]+[(round(.15+.05*i,2),round(.2+.05*i,2)) for i in range(7)]
    for za,zb in windows:
        base=valid&(z>=za)&(z<zb)
        for interior in [False,True]:
            parent=base&(d['BOUNDARY_MPC']>20) if interior else base
            eligible=[]
            for c in np.unique(cells[parent]):
                ok=True
                for sel in [parent&(cells==c),parent&low&(cells==c)]:
                    for k in range(4):
                        w=p[sel,k];eff=w.sum()**2/(w@w) if (w@w)>0 else 0
                        ok &= eff>=20
                if ok:eligible.append(int(c))
            common=parent&np.isin(cells,eligible)
            record=dict(z=[za,zb],boundary_over_20_mpc=interior,parent_rows=int(parent.sum()),common_parent_rows=int(common.sum()),common_low_rows=int((common&low).sum()),cells=eligible,arms=[])
            if eligible:
                for selected in [False,True]:
                    sel=common&low if selected else common; values=[];bs=[];weights=[]
                    for c in eligible:
                        ii=np.flatnonzero(sel&(cells==c)); vals=[];bb=[]
                        for k in range(4):
                            den=np.bincount(blocks[ii],weights=p[ii,k],minlength=nb);num=np.bincount(blocks[ii],weights=p[ii,k]*q[ii],minlength=nb); vals.append(num.sum()/den.sum());db=bw@den;bb.append(np.divide(bw@num,db,out=np.full(200,np.nan),where=db>0))
                        values.append(vals);bs.append(np.array(bb).T);weights.append((common&(cells==c)).sum())
                    w=np.array(weights)/sum(weights);v=w@np.array(values); b=np.einsum('i,ijk->jk',w,np.array(bs));diff=v[3]-v[0];ci=np.nanquantile(b[:,3]-b[:,0],[.025,.975])
                    record['arms'].append(dict(low_width=selected,class_quenched=v.tolist(),contrast_pp=float(diff*100),ci95_pp=(ci*100).tolist()))
            # Cell-fixed-effect partial Spearman: ranks globally within arm, center each cell,
            # residualize continuous ranked mass,z within cells; no class-support requirement.
            record['selection_profiles']=[]
            for selected in [False,True]:
                ss=parent&low if selected else parent
                record['selection_profiles'].append(dict(low_width=selected,n=int(ss.sum()),mean_pweb=p[ss].mean(0).tolist(),trace_quantiles=np.quantile(eig[ss].sum(1),[.1,.5,.9]).tolist(),width_quantiles=np.quantile(width[ss],[.1,.5,.9]).tolist()))
            record['partial_rank']=[]
            for selected in [False,True]:
                sel=parent&low if selected else parent
                counts={c:int((sel&(cells==c)).sum()) for c in np.unique(cells[sel])}
                ii=np.flatnonzero(sel&np.isin(cells,[c for c,n in counts.items() if n>=40]))
                if len(ii)<100:continue
                a=np.column_stack([eig[ii,0],eig[ii].sum(1),q[ii],np.log10(np.maximum(sfr[ii]/m[ii],1e-14)),mass[ii],z[ii]])
                a=np.column_stack([rankdata(a[:,j]) for j in range(a.shape[1])])
                for c in np.unique(cells[ii]): a[cells[ii]==c]-=a[cells[ii]==c].mean(0)
                control=a[:,4:]; residual=a[:,:4]-control@np.linalg.lstsq(control,a[:,:4],rcond=None)[0]
                corr=np.corrcoef(residual.T)
                record['partial_rank'].append(dict(low_width=selected,n=len(ii),lambda1_quenched=float(corr[0,2]),trace_quenched=float(corr[1,2]),lambda1_logssfr=float(corr[0,3]),trace_logssfr=float(corr[1,3])))
            result['windows'].append(record)
    dump('test3.json',result)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5)); rows=[r for r in result['windows'] if not r['boundary_over_20_mpc'] and r['z'][1]-r['z'][0]<.051]
    for selected,color,label in [(False,'#167d9a','All eligible'),(True,'#d66e27','Lowest 50% width per cell')]:
        for ax,field in [(axes[0],'contrast'),(axes[1],'partial')]:
            xx=[];yy=[];lo=[];hi=[]
            for r in rows:
                arms=r['arms'] if field=='contrast' else r['partial_rank'];arm=next((a for a in arms if a['low_width']==selected),None)
                if arm is None:continue
                xx.append(np.mean(r['z']));yy.append(arm['contrast_pp'] if field=='contrast' else arm['trace_logssfr'])
                if field=='contrast':lo.append(yy[-1]-arm['ci95_pp'][0]);hi.append(arm['ci95_pp'][1]-yy[-1])
            if field=='contrast':ax.errorbar(xx,yy,yerr=[lo,hi],fmt='o-',capsize=3,color=color,label=label)
            else:ax.plot(xx,yy,'o-',color=color,label=label)
    axes[0].set_ylabel('Soft cluster − void quenched contrast (pp)');axes[1].set_ylabel('Partial rank correlation: trace vs log sSFR')
    for ax in axes:ax.set_xlabel('Redshift');ax.axhline(0,color='gray',lw=.5);ax.legend(fontsize=8)
    fig.suptitle('Test 3: 0.2 dex × 0.05 redshift controls; actual eigenvalue interval width');fig.tight_layout();fig.savefig(OUT/'test3_uncertainty.png',dpi=160);plt.close(fig)
    dump('VALIDATION.json',dict(unique_ids=True,probability_normalization_checked=True,finite_outcomes_checked=True,uncertainty_half_selection_nested=True,common_class_support_minimum_kish=20,bootstrap_draws=200,bootstrap_sky_blocks=int(nb),scope='conditional spatial sampling uncertainty only; rank correlations descriptive, no calibration claim'))
    print('Test 3 complete',flush=True)
if __name__=='__main__':
    test1();test3()
