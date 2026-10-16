"""Colour-split Loa properties; run only in an approved CPU allocation."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[k]='1'
import json
import socket
from pathlib import Path
import numpy as np
from run import OUT as SOURCE, PHOT, controlled, write_json, sha
OUT=SOURCE.parent/'loa_colour_split_20261016_v1'
LABELS={'quenched':'Quenched fraction (sSFR < 10⁻¹¹ yr⁻¹)', 'logmass':'Mean log₁₀ stellar mass [M☉]',
        'logssfr':'Mean log₁₀ sSFR [yr⁻¹; floor −14]', 'logsfr':'Mean log₁₀ SFR [M☉/yr; floor −14]',
        'colour':'Observed g − r [mag]', 'logage':'Mean log₁₀ fitted age [Gyr]',
        'av':'Fitted dust Aᵥ [mag]', 'fagn':'Fitted AGN contribution'}
def masks(c):
    return {'red':np.isfinite(c)&(c>.75)&(c<1.5), 'blue':np.isfinite(c)&(c<.75)}
def threshold(curve,hi):
    bins=[b for b in curve if b['zlo']>=.2 and b['zhi']<=hi]
    assert len(bins)==round((hi-.2)*100) and all(b['usable'] for b in bins)
    return float(np.ceil(100*max(b['logmass_limit'] for b in bins))/100)
def colour_limits():
    """Independent colour-population completeness audit from the Bright parent."""
    import fitsio
    b=np.load(SOURCE/'bright_parent.npz'); order=np.argsort(b['TARGETID']); ids=b['TARGETID'][order]
    col=np.full(len(ids),np.nan); seen=np.zeros(len(ids),bool)
    with fitsio.FITS(PHOT) as f:
        for start in range(0,f[1].get_nrows(),300000):
            a=f[1].read(rows=np.arange(start,min(start+300000,f[1].get_nrows())),columns=['TARGETID','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R','ZWARN','DELTACHI2','Z_not4clus'])
            at=np.searchsorted(ids,a['TARGETID']); ok=(at<len(ids))&(ids[np.minimum(at,len(ids)-1)]==a['TARGETID'])
            ok&=(a['ZWARN']==0)&(a['DELTACHI2']>40)&(a['Z_not4clus']>0)&(a['Z_not4clus']<.55)
            ix=order[at[ok]]; assert not seen[ix].any();seen[ix]=True
            x=a[ok]
            with np.errstate(divide='ignore',invalid='ignore'):
                c=-2.5*np.log10((x['FLUX_G']/x['MW_TRANSMISSION_G'])/(x['FLUX_R']/x['MW_TRANSMISSION_R']))
            valid=(x['FLUX_G']>0)&(x['FLUX_R']>0)&(x['MW_TRANSMISSION_G']>0)&(x['MW_TRANSMISSION_R']>0)
            col[ix]=np.where(valid,c,np.nan)
    assert seen.all()
    curve=json.loads((SOURCE/'mass_limits.json').read_text()); groups=masks(col); audit=[]
    offset=b['r']-b['rlim']
    for old in curve:
        if old['zlo']<.2 or old['zhi']>.3:continue
        rec=dict(old);rec['colour_populations']={};limits=[old['logmass_limit']]
        for group,mask in groups.items():
            ii=np.flatnonzero(mask&(b['z']>=old['zlo'])&(b['z']<old['zhi']))
            item={'rows':len(ii),'usable':len(ii)>=100}
            if item['usable']:
                faint=ii[offset[ii]>=np.quantile(offset[ii],.8)]
                values=b['mass'][faint]+.4*offset[faint]
                lim=float(np.quantile(values,.95)); limits.append(lim)
                # Independent sorted order-statistic reconstruction.
                v=np.sort(values); pos=.95*(len(v)-1); j=int(pos)
                check=v[j]+(pos-j)*(v[min(j+1,len(v)-1)]-v[j]);assert abs(check-lim)<1e-10
                item.update(p95=lim,faint_rows=len(faint))
            rec['colour_populations'][group]=item
        rec['usable']=old['usable'] and all(x['usable'] for x in rec['colour_populations'].values())
        rec['logmass_limit']=max(limits) if rec['usable'] else None;audit.append(rec)
    return audit

def analyse(a,mask,width=.1):
    import healpy as hp
    a=a[mask]; m=np.log10(a['MASS_CG_15']);z=a['Z'];p=a['P_WEB'].astype(float)
    assert np.isfinite(p).all() and (p>=0).all() and np.allclose(p.sum(1),1,atol=1e-5)
    groups=masks(a['GR_OBS']);hard=np.eye(4)[p.argmax(1)]
    weights={g+'_'+method:w*keep[:,None] for g,keep in groups.items() for method,w in [('soft',p),('hard',hard)]}
    sf=a['SFR_CG_15']; valid=np.isfinite(sf)&(sf>=0)
    with np.errstate(divide='ignore',invalid='ignore'):
        ss=np.where(valid,np.log10(sf/a['MASS_CG_15']),np.nan)
        props={'quenched':np.where(valid,(ss<-11).astype(float),np.nan),'logmass':m,'logssfr':np.maximum(ss,-14),
               'logsfr':np.where(valid,np.maximum(np.log10(sf),-14),np.nan),'colour':a['GR_OBS'],
               'logage':np.log10(np.where(a['AGE_CG_15']>0,a['AGE_CG_15']/1000,np.nan)),
               'av':np.where(a['AV_CG_15']>=0,a['AV_CG_15'],np.nan),
               'fagn':np.where((a['FAGN_CG_15']>=0)&(a['FAGN_CG_15']<=1),a['FAGN_CG_15'],np.nan)}
    region=hp.ang2pix(8,a['RA'],a['DEC'],lonlat=True);results={}
    for prop,y in props.items():
        use=np.isfinite(y)&(groups['red']|groups['blue']);ii=np.flatnonzero(use)
        zb=np.floor((z[ii]-.2+1e-9)/.01).astype(int)
        cells=zb if prop=='logmass' else np.floor((m[ii]-6)/width).astype(int)*100+zb
        w={k:v[ii] for k,v in weights.items()};r=controlled(y[ii],w,cells,region[ii])
        if r['status']=='measured':
            kept=np.isin(cells,r['common_cell_ids']);r['colour_rows']={g:int(keep[ii][kept].sum()) for g,keep in groups.items()}
            r['matched_z_range']=[float(z[ii][kept].min()),float(z[ii][kept].max())]
            # Independent direct means using recorded common strata and mix.
            for method,ww in w.items():
                check=np.zeros(4)
                for cell,mix in zip(r['common_cell_ids'],r['reference_weights']):
                    q=cells==cell;check+=mix*(ww[q].T@y[ii][q])/ww[q].sum(0)
                assert np.allclose(check,r['methods'][method]['mean'],atol=1e-11)
        results[prop]=r
    return results

def render(results, windows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for name,blocks in results.items():
        fig,axs=plt.subplots(2,4,figsize=(17,8)); block=blocks['0.1']
        for ax,(prop,label) in zip(axs.flat,LABELS.items()):
            rec=block[prop]
            if rec['status']=='measured':
                for group,col in [('red','#b52c3c'),('blue','#176bad')]:
                    for method,style in [('soft','-'),('hard','--')]:
                        v=rec['methods'][group+'_'+method];ci=np.array(v['ci95'])
                        ax.plot(range(4),v['mean'],marker='o' if method=='soft' else 's',ls=style,c=col,label=group+' · '+method)
                        ax.vlines(range(4),ci[0],ci[1],color=col,alpha=.65)
                n=rec['colour_rows'];note=f"Matched red {n['red']:,}; blue {n['blue']:,}"
            else:note='Insufficient shared support';ax.text(.05,.5,note,transform=ax.transAxes)
            if prop=='logmass':note+=' · z control only'
            ax.set_title(label+'\n'+note,fontsize=9);ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'],rotation=20);ax.grid(alpha=.2)
        axs[0,0].legend(fontsize=8)
        w=windows[name];fig.suptitle(f"Loa 0.20 ≤ z < {w['zhi']:.2f} · log₁₀ M★ ≥ {w['logmass']:.2f}\nRed: 0.75 < g−r < 1.5; blue: g−r < 0.75 · shared mass/z support · 95% sky-block intervals")
        fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(OUT/f'properties_{name}.png',dpi=180);fig.savefig(OUT/f'properties_{name}.pdf');plt.close(fig)

def main():
    import fitsio
    from numpy.lib.recfunctions import append_fields
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('Compute allocation required')
    OUT.mkdir(exist_ok=False)
    curve=colour_limits();write_json(OUT/'colour_mass_limits.json',curve)
    a=np.load(SOURCE/'analysis_input.npy',allow_pickle=False);m=np.log10(a['MASS_CG_15']);groups=masks(a['GR_OBS']); p=a['P_WEB'].astype(float)
    windows={};cuts={}
    for hi in [.30,.25,.23,.22]:
        key=f'z0200{round(hi*100):02d}';lim=threshold(curve,hi);cut=(a['Z']>=.2)&(a['Z']<hi)&(m>=lim)
        counts={g:int((cut&v).sum()) for g,v in groups.items()}; hard={g:np.bincount(p[cut&v].argmax(1),minlength=4).tolist() for g,v in groups.items()}
        enough=min(counts.values())>=20000 and min(min(v) for v in hard.values())>=100
        windows[key]={'zhi':hi,'logmass':lim,'mass_msun':float(10**lim),'rows':int(cut.sum()),'counts':counts,'hard_counts':hard,'adequate_counts':enough};cuts[key]=cut
    adequate=[k for k,v in windows.items() if v['zhi']<.3 and v['adequate_counts']]
    chosen=min(adequate,key=lambda k:windows[k]['zhi']) if adequate else None
    write_json(OUT/'window_counts.json',{'windows':windows,'chosen':chosen,'rule':'Smallest tested zmax >=.22 with >=20000 each colour and >=100 each hard class per colour; chosen without property outcomes'})
    results={}
    for key in dict.fromkeys(['z020030']+([chosen] if chosen else [])):
        results[key]={str(width):analyse(a,cuts[key],width) for width in [.1,.05]}
        sub=a[cuts[key]];flag=np.full(len(sub),-1,dtype='i2');cm=masks(sub['GR_OBS']);flag[cm['blue']]=0;flag[cm['red']]=1
        sub=append_fields(sub,'COLOUR_CLASS',flag,usemask=False)
        fitsio.write(OUT/f'LOA_COLOUR_MASS_COMPLETE_{key}.fits',sub,header={'MCFIX':windows[key]['logmass'],'MCZLO':.2,'MCZHI':windows[key]['zhi'],'COLRED':'0.75 < g-r < 1.5','COLBLUE':'g-r < 0.75'})
        fresh=fitsio.read(OUT/f'LOA_COLOUR_MASS_COMPLETE_{key}.fits',rows=np.linspace(0,len(sub)-1,min(257,len(sub)),dtype=int))
        for col in a.dtype.names:assert np.array_equal(fresh[col],sub[col][np.linspace(0,len(sub)-1,min(257,len(sub)),dtype=int)],equal_nan=True) if fresh[col].dtype.kind in 'fc' else np.array_equal(fresh[col],sub[col][np.linspace(0,len(sub)-1,min(257,len(sub)),dtype=int)])
    write_json(OUT/'relationships.json',results);render(results,windows)
    write_json(OUT/'MANIFEST.json',{'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'source_manifest_sha256':sha(SOURCE/'MANIFEST.json'),'code':{str(p):sha(p) for p in [Path(__file__),Path(__file__).with_name('run.py')]},'outputs':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()},'validation':'Independent percentile reconstruction, every reported common-stratum curve, and FITS sample roundtrip passed','limitations':'Empirical flux completeness; observed dereddened colour, no K correction. Conditional sky-bootstrap excludes SED, selection and environment-transfer systematics. No causal inference.'})
    print(json.dumps({'windows':windows,'chosen':chosen}),flush=True)
if __name__=='__main__':main()
