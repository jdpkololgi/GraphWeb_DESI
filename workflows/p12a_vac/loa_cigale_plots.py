"""Reproduce CIGALE wedge figures and extend them to the full matched Loa survey."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
import fitsio,healpy as hp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
import property_environment as old
from property_environment import estimate,NAMES,COLORS,digest

def main(a):
    a.out.mkdir(parents=True,exist_ok=False)
    d=fitsio.read(a.join);d=d.astype(d.dtype.newbyteorder('='));p=d['P_WEB'].astype(np.float64);z=d['Z'];supported=d['SUPPORTED']
    matched=d['CIGALE_MATCHED'];gal=np.char.strip(d['CIGALE_SPECTYPE'].astype('U'))=='GALAXY'
    base=matched&supported&gal
    properties={};valid={}
    for fit in ['15','5']:
        m=d[f'MASS_CG_{fit}'];s=d[f'SFR_CG_{fit}'];ok=base&np.isfinite(m)&np.isfinite(s)&(m>0)&(s>0)
        lm=np.full(len(d),np.nan);ss=lm.copy();lm[ok]=np.log10(m[ok]);ss[ok]=np.log10(s[ok])-lm[ok]
        properties[fit]=(lm,ss);valid[fit]=ok
    region=hp.ang2pix(8,d['RA'],d['DEC'],lonlat=True)
    zbin=np.clip(np.floor((z-.15)/.025).astype(int),0,15)
    wedge=(d['RA']>=120)&(d['RA']<160)&(d['DEC']>=14.5)&(d['DEC']<30.6)&(z>=.2)&(z<.3)
    low=(z>=.2)&(z<.3)
    # This cache contains property columns only, never historical environment predictions.
    kk=matched&gal&wedge
    cache=a.out/'matched_wedge_properties.parquet'
    pd.DataFrame({'TARGETID':d['TARGETID'][kk],'ra':d['RA'][kk],'dec':d['DEC'][kk],'z':z[kk],
                  'MASS_CG':d['MASS_CG_15'][kk],'SFR_CG':d['SFR_CG_15'][kk],
                  'MASS_CG5':d['MASS_CG_5'][kk],'SFR_CG5':d['SFR_CG_5'][kk]}).to_parquet(cache,index=False)
    old.CIG=cache;old.main(a.out/'original_wedge')
    wr=a.out/'original_wedge'/'RESULTS.json'
    wj=json.loads(wr.read_text());wj['selection']='Unique coordinate/redshift-consistent CIGALE match, SPECTYPE GALAXY, supported VAC, wedge, positive finite mass/SFR. No legacy goodPhoto cut or new SED-quality cut.'
    wr.write_text(json.dumps(wj,indent=2)+'\n')
    wm=a.out/'original_wedge'/'ARTIFACT_SHA256.json'
    wm.write_text(json.dumps({f.name:digest(f) for f in wm.parent.iterdir() if f.is_file() and f!=wm},indent=2)+'\n')
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
    footer='Provisional Loa environments; real-DESI coverage unverified. High-z population mismatch; trends are not physical evolution.\nCIGALE SFRs may be affected by Tractor modelling, especially at low z. Descriptive probability-weighted relationships.'
    pdf=PdfPages(a.out/'FULL_SURVEY_CIGALE_environment.pdf');files=[]
    def save(fig,name):
        fig.text(.5,.025,footer,ha='center',va='bottom',fontsize=8,color='#555555')
        fig.savefig(a.out/f'{name}.png',dpi=180);pdf.savefig(fig);plt.close(fig);files.append(f'{name}.png')
    distributions={}
    for fit,scope,mask in [('15','full_range',np.ones(len(d),bool)),('15','lowz_full_sky',low),('5','lowz_full_sky',low)]:
        k=valid[fit]&mask;lm,ss=properties[fit]
        fig,axs=plt.subplots(2,2,figsize=(12,9));fig.subplots_adjust(left=.08,right=.86,bottom=.17,top=.85,wspace=.22,hspace=.25)
        fig.suptitle(f'Stellar mass and sSFR · CIGALE {fit}-band · {scope.replace("_"," ")}',fontsize=17)
        fig.text(.5,.90,f'{int(k.sum()):,} galaxies; same sample in each panel, weighted by current environment probabilities\nEach class is normalized by its total probability weight; shared density scale.',ha='center',fontsize=10)
        xe=np.linspace(8.5,12.5,65);ye=np.linspace(-14,-8,81);report=[]
        for j,ax in enumerate(axs.flat):
            h,_,_=np.histogram2d(lm[k],ss[k],bins=[xe,ye],weights=p[k,j]);den=p[k,j].sum();rho=h/(den*np.diff(xe)[:,None]*np.diff(ye)[None,:])
            im=ax.pcolormesh(xe,ye,rho.T,norm=LogNorm(vmin=.001,vmax=2),cmap='magma',rasterized=True)
            ax.axhline(-11,color='#44c6ce',ls='--',lw=1);ax.set(title=NAMES[j],xlabel='log₁₀ stellar mass [M☉]',ylabel='log₁₀ sSFR [yr⁻¹]')
            report.append({'class':NAMES[j],'weight':float(den),'fraction_in_axes':float(h.sum()/den)})
        fig.colorbar(im,cax=fig.add_axes([.89,.23,.018,.5]),label='Probability-weighted density [dex⁻²]')
        name=f'distribution_cg{fit}_{scope}';distributions[name]=report;save(fig,name)
    # Use identical galaxies, cell definitions and reference weights for the two property fits.
    m15,ss15=properties['15'];m5,ss5=properties['5']
    common=valid['15']&valid['5']&(m15>=8)&(m15<13)
    masscell=np.zeros(len(d),dtype=int);masscell[common]=np.floor((m15[common]-8)/.1).astype(int)
    def est(y,mask,mc=True):
        k=common&mask&np.isfinite(y);cells=zbin[k]+16*masscell[k] if mc else zbin[k]
        _,cells=np.unique(cells,return_inverse=True);_,reg=np.unique(region[k],return_inverse=True)
        return estimate(y[k],p[k],cells,reg,min_eff=20,nboot=64)
    summaries={}
    for scope,mask in [('full_range',np.ones(len(d),bool)),('lowz_full_sky',low)]:
        fig,axs=plt.subplots(1,3,figsize=(15,5.5));fig.subplots_adjust(left=.07,right=.98,bottom=.25,top=.76,wspace=.3)
        fig.suptitle('Galaxy properties versus environment · '+scope.replace('_',' '),fontsize=18)
        fig.text(.5,.84,'CG15 and CG5 use identical galaxies and control cells. Mass: common redshift mix.\nsSFR / low-sSFR: common CG15-mass and redshift mix. Error bars: 16–84% sky-block bootstrap.',ha='center',fontsize=10)
        summaries[scope]={}
        for j,(metric,title,mc) in enumerate([('mass','Mean log₁₀ stellar mass [M☉]',False),('logssfr','Mean log₁₀ sSFR [yr⁻¹]',True),('low_ssfr','Low-sSFR fraction (< −11 dex)',True)]):
            summaries[scope][metric]={}
            for fit,ls,offset in [('15','-',-.04),('5','--',.04)]:
                lm,ss=properties[fit];y=lm if metric=='mass' else ss if metric=='logssfr' else (ss<-11).astype(float)
                r=est(y,mask,mc);summaries[scope][metric][fit]=r
                if r:
                    x=np.arange(4)+offset;ci=np.array(r['interval']);axs[j].plot(x,r['mean'],ls,label=f'CG{fit}',color='#333333');axs[j].vlines(x,ci[0],ci[1],color=COLORS);axs[j].scatter(x,r['mean'],c=COLORS,marker='o' if fit=='15' else 's')
            axs[j].set(xticks=np.arange(4),xticklabels=NAMES,ylabel=title);axs[j].legend()
        save(fig,f'controlled_{scope}')
    fig,axs=plt.subplots(1,2,figsize=(13,5.5));fig.subplots_adjust(left=.09,right=.98,bottom=.25,top=.80,wspace=.28)
    fig.suptitle('Fixed stellar mass · CG15 · 0.20 ≤ z < 0.30, full sky',fontsize=18)
    fig.text(.5,.86,'Each mass bin shares a redshift distribution across environments; shading: 16–84% spatial bootstrap.',ha='center',fontsize=10)
    bymass={};edges=np.arange(9,12.51,.5);centres=(edges[:-1]+edges[1:])/2
    for ax,(metric,y,label) in zip(axs,[('logssfr',ss15,'Mean log₁₀ sSFR [yr⁻¹]'),('low_ssfr',(ss15<-11).astype(float),'Low-sSFR fraction')]):
        rr=[est(y,low&(m15>=lo)&(m15<hi),False) for lo,hi in zip(edges[:-1],edges[1:])];bymass[metric]=rr
        for j in range(4):
            mu=[r['mean'][j] if r else np.nan for r in rr];lo=[r['interval'][0][j] if r else np.nan for r in rr];hi=[r['interval'][1][j] if r else np.nan for r in rr]
            ax.plot(centres,mu,marker=['o','s','^','D'][j],color=COLORS[j],label=NAMES[j]);ax.fill_between(centres,lo,hi,color=COLORS[j],alpha=.12)
        ax.set(xlabel='CG15 log₁₀ stellar mass [M☉]',ylabel=label);ax.legend(fontsize=9)
    save(fig,'fixed_mass_lowz_full_sky')
    shells={};fig,axs=plt.subplots(1,2,figsize=(13,5.5));fig.subplots_adjust(left=.08,right=.98,bottom=.25,top=.76,wspace=.3)
    fig.suptitle('Environment contrasts by redshift · selection-sensitive diagnostics',fontsize=17)
    fig.text(.5,.84,'Knot minus void, with common CG15-mass / redshift controls. Same galaxies for CG15 and CG5.\nThese are not measurements of physical evolution; model-transfer and CIGALE systematics are not in the intervals.',ha='center',fontsize=10)
    zb=np.array([.15,.25,.35,.45,.55]);centres=(zb[1:]+zb[:-1])/2
    for ax,metric in zip(axs,['logssfr','low_ssfr']):
        shells[metric]={}
        for fit in ['15','5']:
            ss=properties[fit][1];y=ss if metric=='logssfr' else (ss<-11).astype(float)
            rr=[est(y,(z>=lo)&(z<hi)) for lo,hi in zip(zb[:-1],zb[1:])];shells[metric][fit]=rr
            ax.plot(centres,[r['knot_minus_void'] if r else np.nan for r in rr],'o-' if fit=='15' else 's--',label=f'CG{fit}')
            ax.fill_between(centres,[r['contrast_interval'][0] if r else np.nan for r in rr],[r['contrast_interval'][1] if r else np.nan for r in rr],alpha=.15)
        ax.axhline(0,color='gray',ls=':');ax.set(xlabel='Redshift',ylabel='Knot − void '+('mean log sSFR [dex]' if metric=='logssfr' else 'low-sSFR fraction'));ax.legend()
    save(fig,'redshift_contrast_sensitivity')
    fig,axs=plt.subplots(1,2,figsize=(13,5.3));fig.subplots_adjust(left=.08,right=.98,bottom=.25,top=.79,wspace=.28)
    fig.suptitle('CIGALE matching and property availability',fontsize=18)
    fig.text(.5,.85,'Denominator: all supported VAC galaxies in each redshift bin, weighted by environment probability.',ha='center',fontsize=10)
    edges=np.arange(.15,.551,.05);xc=(edges[:-1]+edges[1:])/2;completeness={'redshift_edges':edges.tolist(),'series':{}}
    for ax,mask,title in zip(axs,[matched,valid['15']],['Unique coordinate-consistent match','Usable CG15 galaxy mass and SFR']):
        completeness['series'][title]={}
        for j in range(4):
            rates=[]
            for lo,hi in zip(edges[:-1],edges[1:]):
                k=supported&(z>=lo)&(z<hi);rates.append(float(p[k&mask,j].sum()/p[k,j].sum()))
            completeness['series'][title][NAMES[j]]=rates
            ax.plot(xc,rates,marker=['o','s','^','D'][j],color=COLORS[j],label=NAMES[j])
        ax.set(xlabel='Redshift',ylabel='Fraction of supported VAC',ylim=(0,1.02),title=title);ax.legend(fontsize=9)
    save(fig,'matching_completeness');pdf.close()
    availability={}
    for label,mask in [('all',np.ones(len(d),bool)),('lowz_full_sky',low),('original_wedge',wedge)]:
        den=p[supported&mask].sum(axis=0);availability[label]={'supported_rows':int((supported&mask).sum()),'matched_rows':int((matched&mask).sum()),'valid15':int((valid['15']&mask).sum()),'valid5':int((valid['5']&mask).sum()),'cg15_class_availability':(p[valid['15']&mask].sum(axis=0)/den).tolist()}
    report={'join':str(a.join),'join_sha256':digest(a.join),'script_sha256':digest(__file__),
            'matched_spectypes':dict(zip(*[x.tolist() for x in np.unique(d['CIGALE_SPECTYPE'][matched].astype('U'),return_counts=True)])),
            'availability':availability,'completeness_by_redshift':completeness,'controlled':summaries,'fixed_mass':bymass,'shells':shells,'distributions':distributions,'figures':files,
            'analysis_selection':'Unique position/redshift-consistent CIGALE match, supported VAC, CIGALE SPECTYPE GALAXY, positive finite mass/SFR. Controlled fit comparison further requires both fits valid and 8<=CG15 logM<13. No new VAC cuts; no SED-quality flag supplied; CHI2 is not used as CIGALE quality.',
            'controls':'CG15 mass bins 0.1 dex; z bins .025; >=20 effective rows in every class; fixed common reference;64 nside8 bootstrap draws. CG5 is conditioned on CG15 mass for identical comparison support. No property error covariance or model transfer systematics propagated.'}
    (a.out/'RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'availability':availability,'completeness_by_redshift':completeness,'controlled':summaries},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--join',type=Path,required=True);p.add_argument('--out',type=Path,required=True);main(p.parse_args())
