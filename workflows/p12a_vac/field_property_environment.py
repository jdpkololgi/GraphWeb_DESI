"""Matched full-Loa CIGALE relations: frozen P12-A versus local CFM marginals.

No original VAC mutation. Independently sampled atlas tiles are not a global
joint posterior; bootstrap intervals below are selected-sample diagnostics.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages

NAMES=['Void','Sheet','Filament','Knot']
CLASS_COLORS=['#2585bb','#aa8310','#d25b42','#8954af']
MODEL_NAMES=['P12-A','CFM field marginals']
MODEL_COLORS=['#2864a0','#ce7219']
FOOTER=('Exploratory, selected-sample associations; neither model is calibrated on real Loa. '
        'Known high-z transfer and CIGALE-selection limitations.\n'
        'Atlas tiles are independent local conditionals, not a joint full-survey posterior. '
        'Intervals omit SED and model-transfer systematics.')


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def native(v):
    if isinstance(v,dict):return {str(k):native(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)):return [native(x) for x in v]
    if isinstance(v,np.ndarray):return native(v.tolist())
    if isinstance(v,np.generic):return native(v.item())
    if isinstance(v,float) and not np.isfinite(v):return None
    return v


def paired_estimate(y,weights,cells,regions,min_eff=20,nboot=64):
    """One retained cell set and reference mixture for every compared model."""
    if not len(y):return None
    values,ci=np.unique(cells,return_inverse=True);_,ri=np.unique(regions,return_inverse=True)
    nc,nr=len(values),int(ri.max())+1;nm=len(weights)
    count=np.zeros((nm,4,nc,nr));sums=np.zeros_like(count);squared=np.zeros_like(count)
    flat=ci*nr+ri
    for m,p in enumerate(weights):
        if p.shape!=(len(y),4) or np.any(p<0) or not np.allclose(p.sum(1),1,atol=2e-6):
            raise ValueError('invalid probabilities for matched estimator')
        for k in range(4):
            count[m,k]=np.bincount(flat,weights=p[:,k],minlength=nc*nr).reshape(nc,nr)
            sums[m,k]=np.bincount(flat,weights=p[:,k]*y,minlength=nc*nr).reshape(nc,nr)
            squared[m,k]=np.bincount(flat,weights=p[:,k]**2,minlength=nc*nr).reshape(nc,nr)
    n=count.sum(-1);eff=n*n/np.maximum(squared.sum(-1),1e-30)
    common=np.all(eff>=min_eff,axis=(0,1))
    if not common.any():return None
    keep=common[ci];ref=np.bincount(ci,minlength=nc)[common].astype(float);ref/=ref.sum()
    def calc(c,t):return np.sum(t[:,:,common]/c[:,:,common]*ref,axis=-1)
    def raw(c,t):return t[:,:,common].sum(-1)/c[:,:,common].sum(-1)
    mu=calc(n,sums.sum(-1));rm=raw(n,sums.sum(-1))
    rng=np.random.default_rng(20260930);boots=[];rawboots=[]
    if nr>=2:
        for _ in range(nboot):
            mult=np.bincount(rng.integers(0,nr,nr),minlength=nr)
            nn=count@mult;ss=sums@mult
            if np.any(nn[:,:,common]<=0):continue
            boots.append(calc(nn,ss));rawboots.append(raw(nn,ss))
    ok=len(boots)>=max(2,.8*nboot)
    interval=np.quantile(boots,[.16,.84],axis=0) if ok else np.full((2,nm,4),np.nan)
    raw_interval=np.quantile(rawboots,[.16,.84],axis=0) if ok else interval.copy()
    contrast=mu[:,3]-mu[:,0]
    cb=np.asarray(boots)[:,:,3]-np.asarray(boots)[:,:,0] if ok else None
    return dict(mean=mu,interval=interval,raw_mean=rm,raw_interval=raw_interval,
        contrast=contrast,contrast_interval=np.quantile(cb,[.16,.84],axis=0) if ok else np.full((2,nm),np.nan),
        paired_contrast_difference=float(contrast[1]-contrast[0]) if nm==2 else None,
        paired_difference_interval=np.quantile(cb[:,1]-cb[:,0],[.16,.84]) if ok and nm==2 else None,
        rows=len(y),common_rows=int(keep.sum()),common_cells=int(common.sum()),sky_blocks=nr,
        bootstrap_replicates=len(boots),bootstrap_interval_available=bool(ok),
        cell_ids=values[common],reference_weights=ref,
        effective_rows=np.array([p[keep].sum(0)**2/np.maximum((p[keep]**2).sum(0),1e-30) for p in weights]))


def fixed_reference(y,p,cells,estimate):
    """Monte Carlo sensitivity on frozen cells/reference, no selection re-fitting."""
    ids=np.asarray(estimate['cell_ids']);ref=np.asarray(estimate['reference_weights'])
    mapping={int(v):i for i,v in enumerate(ids)}
    ix=np.array([mapping.get(int(v),-1) for v in cells]);keep=ix>=0
    count=np.stack([np.bincount(ix[keep],weights=p[keep,j],minlength=len(ids)) for j in range(4)])
    sums=np.stack([np.bincount(ix[keep],weights=p[keep,j]*y[keep],minlength=len(ids)) for j in range(4)])
    if np.any(count==0):return None
    return ((sums/count)*ref).sum(1)


def tests():
    cells=np.repeat([0,1],400);regions=np.arange(800)%13;y=cells*3.
    p=np.vstack([np.tile([.4,.3,.2,.1],(400,1)),np.tile([.1,.2,.3,.4],(400,1))])
    r=paired_estimate(y,[p,p[:,::-1]],cells,regions,min_eff=10,nboot=16)
    assert np.allclose(r['mean'],1.5) and not np.allclose(r['raw_mean'][0],r['raw_mean'][1])
    r=paired_estimate(np.ones(800)*7,[p,p[:,::-1]],cells,regions,min_eff=10,nboot=16)
    assert np.allclose(r['mean'],7) and np.allclose(r['raw_mean'],7)
    hard=np.tile(np.eye(4),(200,1));y=np.tile(np.arange(4),200)+cells*3
    r=paired_estimate(y,[hard,hard],cells,regions,min_eff=10,nboot=16)
    assert np.allclose(r['mean'],np.arange(4)+1.5) and r['paired_contrast_difference']==0
    np.testing.assert_allclose(fixed_reference(y,hard,cells,r),r['mean'][0])
    absent=np.tile([.5,.5,0,0],(800,1))
    assert paired_estimate(y,[p,absent],cells,regions,min_eff=10,nboot=16) is None
    print('Four paired-estimator/control/support checks passed',flush=True)


def main(a):
    if 'SLURM_JOB_ID' not in os.environ:raise RuntimeError('catalogue analysis requires a Slurm allocation')
    a.out.mkdir(parents=True,exist_ok=False)
    plan=json.loads((a.root/'PLAN.json').read_text());atlas=json.loads((a.root/'ATLAS_COMPLETE.json').read_text())
    if not atlas['passed'] or digest(a.root/'vac_properties.npz')!=plan['properties_sha256']:
        raise ValueError('atlas/property provenance mismatch')
    if digest(atlas['output']['path'])!=atlas['output']['sha256']:raise ValueError('field marginals changed')
    with np.load(a.root/'vac_properties.npz') as f:d={k:f[k] for k in f.files}
    with np.load(atlas['output']['path']) as f:field={k:f[k] for k in f.files}
    if len(d['TARGETID'])!=atlas['rows'] or len(np.unique(d['TARGETID']))!=atlas['rows']:
        raise ValueError('VAC row identity not unique')
    z=d['Z'];p12=d['P_WEB'].astype(float);pf=field['p_web'].astype(float)
    gal=np.char.strip(d['CIGALE_SPECTYPE'].astype('U'))=='GALAXY'
    base=d['SUPPORTED']&d['CIGALE_MATCHED']&gal;props={};valid={}
    for fit in ('15','5'):
        mass=d['MASS_CG_'+fit];sfr=d['SFR_CG_'+fit]
        ok=base&np.isfinite(mass)&np.isfinite(sfr)&(mass>0)&(sfr>0)
        lm=np.full(len(z),np.nan);ss=lm.copy();lm[ok]=np.log10(mass[ok]);ss[ok]=np.log10(sfr[ok])-lm[ok]
        props[fit]=(lm,ss);valid[fit]=ok
    common=valid['15']&valid['5']&field['field_eligible']&(props['15'][0]>=8)&(props['15'][0]<13)
    regions=hp.ang2pix(8,d['RA'],d['DEC'],lonlat=True)
    zb=np.clip(np.floor((z-.15)/.025).astype(int),0,15)
    mb=np.zeros(len(z),int);mb[common]=np.floor((props['15'][0][common]-8)/.1).astype(int)
    low=(z>=.2)&(z<.3);weights=[p12,pf]
    def estimate(y,scope,mass_control=True):
        k=common&scope&np.isfinite(y);cells=zb[k]+16*mb[k] if mass_control else zb[k]
        return paired_estimate(y[k],[p[k] for p in weights],cells,regions[k])
    report=dict(rows=len(z),usable_cg15=int(valid['15'].sum()),field_eligible=int(field['field_eligible'].sum()),
        common_both_fits=int(common.sum()),common_lowz=int((common&low).sum()),
        class_threshold=.2,draws=atlas['draws'],models=MODEL_NAMES,source_product=plan['property_product'],
        atlas_receipt_sha256=digest(a.root/'ATLAS_COMPLETE.json'),script_sha256=digest(__file__),
        method='Matched galaxies and common cells/reference across models;64 nside8 sky-block resamples. Not joint-field credible intervals.',
        qualification=FOOTER,scopes={},figures=[],split_half={})
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white'})
    pdf=PdfPages(a.out/'CFM_P12A_CIGALE_comparison.pdf')
    def save(fig,name):
        fig.text(.5,.018,FOOTER,ha='center',va='bottom',fontsize=8,color='#555555')
        fig.savefig(a.out/(name+'.png'),dpi=170);pdf.savefig(fig);plt.close(fig);report['figures'].append(name+'.png')
    metrics=[('mass','Mean log₁₀ stellar mass [M☉]',False),('logssfr','Mean log₁₀ sSFR [yr⁻¹]',True),
             ('low_ssfr','Low-sSFR fraction (< −11)',True)]
    for scope_name,scope in [('lowz_full_footprint',low),('full_selected_range',np.ones(len(z),bool))]:
        report['scopes'][scope_name]={}
        fig,axes=plt.subplots(2,3,figsize=(15,8));fig.subplots_adjust(left=.07,right=.98,top=.84,bottom=.16,hspace=.34,wspace=.30)
        fig.suptitle('Matched CIGALE CG15 properties · '+scope_name.replace('_',' '),fontsize=17)
        fig.text(.5,.89,'Identical galaxies, control cells and reference mixture across models.\nTop: raw on retained cells; bottom: common z mix (mass), common mass + z mix (sSFR/fraction).',ha='center',fontsize=10)
        for j,(metric,label,mc) in enumerate(metrics):
            report['scopes'][scope_name][metric]={}
            for fit in ('15','5'):
                lm,ss=props[fit];y=lm if metric=='mass' else ss if metric=='logssfr' else (ss<-11).astype(float)
                r=estimate(y,scope,mc);report['scopes'][scope_name][metric][fit]=r
                if fit!='15' or r is None:continue
                for i,prefix in enumerate(('raw_','')):
                    for model in range(2):
                        x=np.arange(4)+(-.035 if model==0 else .035);mu=r[prefix+'mean'][model];ci=r[prefix+'interval'][:,model]
                        axes[i,j].plot(x,mu,'o-' if model==0 else 's--',color=MODEL_COLORS[model],label=MODEL_NAMES[model])
                        axes[i,j].vlines(x,ci[0],ci[1],color=MODEL_COLORS[model],lw=1.6)
                    axes[i,j].set(xticks=np.arange(4),xticklabels=NAMES,ylabel=label)
                    if j==0:axes[i,j].legend(fontsize=9)
                axes[1,j].set_title(f'{r["common_rows"]:,} galaxies · {r["common_cells"]} common cells',fontsize=10)
                if scope_name=='lowz_full_footprint':
                    k=common&scope&np.isfinite(y);cells=zb[k]+16*mb[k] if mc else zb[k]
                    lab=field['classes02'][:,k];half=len(lab)//2
                    pleft=np.stack([(lab[:half]==j).mean(0) for j in range(4)],axis=1)
                    pright=np.stack([(lab[half:]==j).mean(0) for j in range(4)],axis=1)
                    left=fixed_reference(y[k],pleft,cells,r);right=fixed_reference(y[k],pright,cells,r)
                    report['split_half'][metric]=dict(first_half=left,second_half=right,
                        contrast_difference=None if left is None or right is None else float((left[3]-left[0])-(right[3]-right[0])),
                        qualification='Frozen full-draw cells/reference. Independent-tile MC sensitivity, not a credible interval.')
        save(fig,'01_controlled_'+scope_name)
    # Reproduce class-weighted bimodal distributions with a common density scale.
    k=common&low;lm,ss=props['15'];xe=np.linspace(8.5,12.5,65);ye=np.linspace(-14,-8,81)
    fig,axes=plt.subplots(2,4,figsize=(17,8));fig.subplots_adjust(left=.055,right=.90,top=.82,bottom=.16,wspace=.22,hspace=.27)
    fig.suptitle(f'Mass–sSFR distributions · {int(k.sum()):,} identical galaxies · 0.20 ≤ z < 0.30',fontsize=16)
    fig.text(.5,.88,'CG15; each panel normalized by its class probability weight; no mass/redshift standardization in this figure.',ha='center',fontsize=10)
    hist=[]
    for m,p in enumerate(weights):
        for j in range(4):
            h,_,_=np.histogram2d(lm[k],ss[k],bins=[xe,ye],weights=p[k,j]);den=float(p[k,j].sum())
            im=axes[m,j].pcolormesh(xe,ye,(h/(den*np.diff(xe)[:,None]*np.diff(ye)[None,:])).T,
                                   cmap='magma',norm=LogNorm(.001,2),rasterized=True)
            axes[m,j].axhline(-11,color='#44c6ce',ls='--',lw=.8)
            axes[m,j].set(title=MODEL_NAMES[m]+' · '+NAMES[j],xlabel='log₁₀ mass [M☉]',ylabel='log₁₀ sSFR [yr⁻¹]')
            hist.append(dict(model=MODEL_NAMES[m],environment=NAMES[j],weight=den,fraction_in_axes=float(h.sum()/den)))
    fig.colorbar(im,cax=fig.add_axes([.93,.24,.012,.48]),label='Class-weighted density [dex⁻²]')
    report['histograms']=hist;save(fig,'02_mass_ssfr_distributions')
    # Fixed-mass panels; both models use the same redshift-balanced cells.
    fig,axes=plt.subplots(2,2,figsize=(13,9));fig.subplots_adjust(left=.09,right=.97,top=.84,bottom=.16,hspace=.3,wspace=.27)
    fig.suptitle('Fixed stellar mass · identical galaxies · 0.20 ≤ z < 0.30',fontsize=17)
    edges=np.arange(9,12.51,.5);centres=(edges[:-1]+edges[1:])/2;report['fixed_mass']={}
    for col,(metric,y,label) in enumerate([('logssfr',ss,'Mean log₁₀ sSFR [yr⁻¹]'),('low_ssfr',(ss<-11).astype(float),'Low-sSFR fraction')]):
        rr=[estimate(y,low&(lm>=lo)&(lm<hi),False) for lo,hi in zip(edges[:-1],edges[1:])];report['fixed_mass'][metric]=rr
        for m in range(2):
            for j in range(4):
                mu=[r['mean'][m,j] if r else np.nan for r in rr]
                lo=[r['interval'][0,m,j] if r else np.nan for r in rr];hi=[r['interval'][1,m,j] if r else np.nan for r in rr]
                axes[m,col].plot(centres,mu,marker=['o','s','^','D'][j],color=CLASS_COLORS[j],label=NAMES[j])
                axes[m,col].fill_between(centres,lo,hi,color=CLASS_COLORS[j],alpha=.12)
            axes[m,col].set(title=MODEL_NAMES[m],xlabel='CG15 log₁₀ mass [M☉]',ylabel=label);axes[m,col].legend(fontsize=8)
        lim=(min(ax.get_ylim()[0] for ax in axes[:,col]),max(ax.get_ylim()[1] for ax in axes[:,col]))
        for ax in axes[:,col]:ax.set_ylim(*lim)
    save(fig,'03_fixed_mass_relations')
    # Selection/support completeness, with the same explicit denominator as the VAC.
    fig,axes=plt.subplots(1,2,figsize=(13,5.5));fig.subplots_adjust(left=.08,right=.98,top=.79,bottom=.24,wspace=.3)
    fig.suptitle('Selection and field-support completeness',fontsize=17)
    fig.text(.5,.85,'Denominator: all supported P12-A VAC galaxies, weighted by their P12-A class probabilities.',ha='center',fontsize=10)
    ze=np.arange(.15,.551,.05);zc=(ze[:-1]+ze[1:])/2;report['completeness']={}
    for ax,name,mask in zip(axes,['Usable CG15 properties','Common field/property analysis sample'],[valid['15'],common]):
        report['completeness'][name]={}
        for j in range(4):
            rates=[]
            for lo,hi in zip(ze[:-1],ze[1:]):
                kk=d['SUPPORTED']&(z>=lo)&(z<hi);rates.append(float(p12[kk&mask,j].sum()/p12[kk,j].sum()))
            report['completeness'][name][NAMES[j]]=rates
            ax.plot(zc,rates,marker=['o','s','^','D'][j],color=CLASS_COLORS[j],label=NAMES[j])
        ax.set(title=name,xlabel='Redshift',ylabel='Retained probability weight',ylim=(0,1.03));ax.legend(fontsize=9)
    save(fig,'04_selection_completeness')
    # Redshift and fit sensitivity are descriptive, never physical evolution.
    fig,axes=plt.subplots(1,2,figsize=(13,5.5));fig.subplots_adjust(left=.08,right=.98,top=.80,bottom=.24,wspace=.28)
    fig.suptitle('Knot − void contrasts by redshift · selection-sensitive',fontsize=16)
    ze=np.array([.15,.25,.35,.45,.55]);zc=(ze[:-1]+ze[1:])/2;report['shells']={}
    for ax,(metric,label) in zip(axes,[('logssfr','Mean log sSFR contrast [dex]'),('low_ssfr','Low-sSFR fraction contrast')]):
        y=ss if metric=='logssfr' else (ss<-11).astype(float)
        rr=[estimate(y,(z>=lo)&(z<hi)) for lo,hi in zip(ze[:-1],ze[1:])];report['shells'][metric]=rr
        for m in range(2):
            ax.plot(zc,[r['contrast'][m] if r else np.nan for r in rr],'o-' if m==0 else 's--',color=MODEL_COLORS[m],label=MODEL_NAMES[m])
            ax.fill_between(zc,[r['contrast_interval'][0,m] if r else np.nan for r in rr],
                             [r['contrast_interval'][1,m] if r else np.nan for r in rr],color=MODEL_COLORS[m],alpha=.12)
        ax.axhline(0,color='gray',ls=':');ax.set(xlabel='Redshift',ylabel=label);ax.legend()
    save(fig,'05_redshift_sensitivity')
    # Same galaxies and CG15 control cells for the two SED fits.
    fig,axes=plt.subplots(1,2,figsize=(13,5.5));fig.subplots_adjust(left=.08,right=.98,top=.80,bottom=.24,wspace=.30)
    fig.suptitle('SED-fit sensitivity · controlled knot − void · 0.20 ≤ z < 0.30',fontsize=16)
    for ax,(metric,label) in zip(axes,[('logssfr','Mean log sSFR contrast [dex]'),('low_ssfr','Low-sSFR fraction contrast')]):
        for m in range(2):
            rr=[report['scopes']['lowz_full_footprint'][metric][fit] for fit in ('15','5')]
            x=np.arange(2)+(m-.5)*.08
            ax.plot(x,[r['contrast'][m] if r else np.nan for r in rr],'o-' if m==0 else 's--',color=MODEL_COLORS[m],label=MODEL_NAMES[m])
            for xx,r in zip(x,rr):
                if r:ax.vlines(xx,*r['contrast_interval'][:,m],color=MODEL_COLORS[m])
        ax.axhline(0,color='gray',ls=':');ax.set(xticks=[0,1],xticklabels=['CG15','CG5'],ylabel=label);ax.legend()
    save(fig,'06_sed_fit_sensitivity')
    # Class-definition sensitivity must not masquerade as a matched P12-A test.
    fig,axes=plt.subplots(1,2,figsize=(13,5.5));fig.subplots_adjust(left=.08,right=.98,top=.80,bottom=.24,wspace=.30)
    fig.suptitle('CFM class-definition sensitivity · 0.20 ≤ z < 0.30',fontsize=16)
    fig.text(.5,.85,'Same galaxies/control cells across the two CFM definitions; P12-A remains defined at λth = 0.2.',ha='center',fontsize=10)
    report['threshold_sensitivity']={}
    for ax,(metric,y,label) in zip(axes,[('logssfr',ss,'Mean log sSFR [yr⁻¹]'),('low_ssfr',(ss<-11).astype(float),'Low-sSFR fraction')]):
        kk=common&low;cells=zb[kk]+16*mb[kk]
        r=paired_estimate(y[kk],[pf[kk],field['p_web_threshold0'][kk]],cells,regions[kk]);report['threshold_sensitivity'][metric]=r
        if r:
            for m,labelth in enumerate(['λth = 0.2','λth = 0']):
                ax.plot(np.arange(4),r['mean'][m],'o-' if m==0 else 's--',color=MODEL_COLORS[m],label=labelth)
                ax.vlines(np.arange(4),r['interval'][0,m],r['interval'][1,m],color=MODEL_COLORS[m])
        ax.set(xticks=np.arange(4),xticklabels=NAMES,ylabel=label);ax.legend()
    save(fig,'07_threshold_sensitivity')
    pdf.close()
    report['cap_rows']={cap:dict(all=int(np.sum(d['CAP']==j)),common=int(np.sum(common&(d['CAP']==j)))) for cap,j in [('NGC',1),('SGC',0)]}
    (a.out/'RESULTS.json').write_text(json.dumps(native(report),indent=2,allow_nan=False)+'\n')
    (a.out/'ARTIFACT_SHA256.json').write_text(json.dumps({p.name:digest(p) for p in a.out.iterdir() if p.is_file()},indent=2)+'\n')
    print(json.dumps(native({k:report[k] for k in ('rows','usable_cg15','field_eligible','common_both_fits','common_lowz','split_half')})),flush=True)


def details(a):
    """Local draw-to-draw variation, not sampling or transfer uncertainty."""
    if 'SLURM_JOB_ID' not in os.environ:raise RuntimeError('Use compute allocation')
    a.out.mkdir(parents=True,exist_ok=False)
    plan=json.loads((a.root/'PLAN.json').read_text())
    if digest(a.root/'vac_properties.npz')!=plan['properties_sha256']:raise ValueError('properties changed')
    with np.load(a.root/'vac_properties.npz') as f:d={k:f[k] for k in f.files}
    fig,axes=plt.subplots(1,2,figsize=(13,6));fig.subplots_adjust(left=.08,right=.98,top=.77,bottom=.24,wspace=.3)
    fig.suptitle('Within-patch conditional draw variation · CG15 log sSFR',fontsize=17)
    fig.text(.5,.85,'Unadjusted local means: dots are draw medians; bars are 16–84% ranges of 16 draws.\nNot mass/z-controlled, not sky-bootstrap intervals; no SED or model-transfer uncertainty.',ha='center',fontsize=10)
    report={}
    for ax,cap in zip(axes,['NGC','SGC']):
        receipt=json.loads((a.root/'details'/('detail_'+cap)/'COMPLETE.json').read_text());item=receipt['outputs'][0]
        if digest(item['path'])!=item['sha256']:raise ValueError('detail fields changed')
        with np.load(item['path']) as f:rows=f['rows'];lab=f['classes02'];support=f['galaxy_support']
        mass=d['MASS_CG_15'][rows];sfr=d['SFR_CG_15'][rows]
        keep=d['SUPPORTED'][rows]&d['CIGALE_MATCHED'][rows]&(np.char.strip(d['CIGALE_SPECTYPE'][rows].astype('U'))=='GALAXY')
        keep&=np.isfinite(mass)&np.isfinite(sfr)&(mass>0)&(sfr>0)&(support>=.5)
        keep&=(d['Z'][rows]>=.20)&(d['Z'][rows]<.30)
        if min(receipt['core_support'])<.25:keep[:]=False
        y=np.log10(sfr[keep])-np.log10(mass[keep]);lab=lab[:,keep];means=np.full((len(lab),4),np.nan);counts=np.zeros_like(means)
        for i,draw in enumerate(lab):
            for j in range(4):
                k=draw==j;counts[i,j]=k.sum()
                if k.sum()>=20:means[i,j]=y[k].mean()
        # Only show a class interval if every draw supplies >=20 galaxies.
        valid=np.all(np.isfinite(means),axis=0);quant=np.full((3,4),np.nan)
        quant[:,valid]=np.quantile(means[:,valid],[.16,.5,.84],axis=0)
        ax.plot(np.arange(4),quant[1],'s',color=MODEL_COLORS[1],label='CFM local draws')
        ax.vlines(np.arange(4),quant[0],quant[2],color=MODEL_COLORS[1],lw=2)
        p=d['P_WEB'][rows][keep];den=p.sum(0);pm=np.sum(p*y[:,None],0)/np.maximum(den,1e-30)
        ax.plot(np.arange(4),pm,'o--',color=MODEL_COLORS[0],label='P12-A weighted mean')
        ax.set(title=f'{cap}: {int(keep.sum()):,} identical owned-core galaxies',xticks=np.arange(4),xticklabels=NAMES,ylabel='Mean log₁₀ sSFR [yr⁻¹]');ax.legend(fontsize=9)
        report[cap]=dict(rows=int(keep.sum()),draws=len(lab),means=means,class_counts=counts,quantile_16_50_84=quant,p12_mean=pm,
                         qualification='Raw local association, all draws require >=20 class members for interval. Not mass/z controlled.')
    lim=(min(ax.get_ylim()[0] for ax in axes),max(ax.get_ylim()[1] for ax in axes))
    for ax in axes:ax.set_ylim(*lim)
    fig.text(.5,.025,FOOTER,ha='center',fontsize=8,color='#555555')
    fig.savefig(a.out/'local_draw_variation.png',dpi=170);fig.savefig(a.out/'local_draw_variation.pdf');plt.close(fig)
    (a.out/'LOCAL_RESULTS.json').write_text(json.dumps(native(report),indent=2,allow_nan=False)+'\n')
    print(json.dumps(native(report)),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path);p.add_argument('--out',type=Path);p.add_argument('--test',action='store_true');p.add_argument('--details',action='store_true')
    a=p.parse_args()
    if a.test:tests()
    else:
        if a.root is None or a.out is None:p.error('--root and --out required')
        details(a) if a.details else main(a)
