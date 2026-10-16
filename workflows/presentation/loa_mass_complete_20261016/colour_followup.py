"""Report sparse-colour support honestly: raw and within-colour controlled curves."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,socket
from pathlib import Path
import numpy as np
from run import OUT as SOURCE,controlled,write_json,sha
from colour_split import OUT as AUDIT,masks,LABELS
OUT=SOURCE.parent/'loa_colour_split_20261016_v2'
def raw(y,w,region):
    _,reg=np.unique(region,return_inverse=True);nr=reg.max()+1
    c=np.array([np.bincount(reg,weights=w[:,j],minlength=nr) for j in range(4)])
    s=np.array([np.bincount(reg,weights=w[:,j]*y,minlength=nr) for j in range(4)])
    mean=s.sum(1)/c.sum(1);boot=[];rng=np.random.default_rng(20261016)
    for _ in range(256):
        mult=np.bincount(rng.integers(nr,size=nr),minlength=nr);den=c@mult
        if (den>0).all():boot.append((s@mult)/den)
    boot=np.array(boot)
    assert np.allclose(mean,w.T@y/w.sum(0))
    return {'mean':mean.tolist(),'ci95':np.quantile(boot,[.025,.975],axis=0).tolist(),'contrast':float(mean[-1]-mean[0]),'contrast_ci95':np.quantile(boot[:,-1]-boot[:,0],[.025,.975]).tolist(),'bootstrap_valid':len(boot)}
def analyse(a):
    import healpy as hp
    results={}
    for group,mask in masks(a['GR_OBS']).items():
        b=a[mask];m=np.log10(b['MASS_CG_15']);z=b['Z'];sf=b['SFR_CG_15'];valid=np.isfinite(sf)&(sf>=0)
        p=b['P_WEB'].astype(float);hard=np.eye(4)[p.argmax(1)];region=hp.ang2pix(8,b['RA'],b['DEC'],lonlat=True)
        with np.errstate(divide='ignore',invalid='ignore'):
            ss=np.where(valid,np.log10(sf/b['MASS_CG_15']),np.nan)
            props={'quenched':np.where(valid,(ss<-11).astype(float),np.nan),'logmass':m,'logssfr':np.maximum(ss,-14),
                   'logsfr':np.where(valid,np.maximum(np.log10(sf),-14),np.nan),'colour':b['GR_OBS'],
                   'logage':np.log10(np.where(b['AGE_CG_15']>0,b['AGE_CG_15']/1000,np.nan)),
                   'av':np.where(b['AV_CG_15']>=0,b['AV_CG_15'],np.nan),
                   'fagn':np.where((b['FAGN_CG_15']>=0)&(b['FAGN_CG_15']<=1),b['FAGN_CG_15'],np.nan)}
        block={}
        for prop,y in props.items():
            ii=np.flatnonzero(np.isfinite(y));rec={'rows':len(ii),'raw':{method:raw(y[ii],w[ii],region[ii]) for method,w in [('soft',p),('hard',hard)]},'controlled':{}}
            for width in [.1,.05]:
                zb=np.floor((z[ii]-.2+1e-9)/.01).astype(int);cells=zb if prop=='logmass' else np.floor((m[ii]-6)/width).astype(int)*100+zb
                # Soft-only support and joint soft/hard support are explicit alternatives.
                for mode,weights in [('soft',{'soft':p[ii]}),('joint',{'soft':p[ii],'hard':hard[ii]})]:
                    r=controlled(y[ii],weights,cells,region[ii])
                    if r['status']=='measured':
                        kept=np.isin(cells,r['common_cell_ids']);r['matched_z_range']=[float(z[ii][kept].min()),float(z[ii][kept].max())]
                        for method,w in weights.items():
                            check=np.zeros(4)
                            for cell,mix in zip(r['common_cell_ids'],r['reference_weights']):
                                use=cells==cell;check+=mix*(w[use].T@y[ii][use])/w[use].sum(0)
                            assert np.allclose(check,r['methods'][method]['mean'],atol=1e-11)
                    rec['controlled'][str(width)+'_'+mode]=r
            block[prop]=rec
        results[group]=block
    return results

def render(results,windows):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    for key,r in results.items():
        fig,axs=plt.subplots(2,4,figsize=(18,8.8))
        for ax,(prop,label) in zip(axs.flat,LABELS.items()):
            notes=[]
            for group,col in [('red','#b52c3c'),('blue','#176bad')]:
                rec=r[group][prop];adj=rec['controlled']['0.1_soft']
                notes.append(f"{group}: {rec['rows']:,} raw / {adj.get('rows',0):,} controlled")
                for method,v,style,marker in [('Raw',rec['raw']['soft'],':','s')]+([('Controlled',adj['methods']['soft'],'-','o')] if adj['status']=='measured' else []):
                    ci=np.array(v['ci95']);mu=np.array(v['mean']);x=np.arange(4)+(-.035 if group=='red' else .035)
                    ax.plot(x,mu,marker=marker,ls=style,c=col,label=group+' · '+method);ax.vlines(x,ci[0],ci[1],color=col,alpha=.65)
            ax.set_title(label+'\n'+'\n'.join(notes),fontsize=9);ax.set_xticks(range(4),['Void','Wall','Filament','Cluster'],rotation=20);ax.grid(alpha=.2)
        axs[0,0].legend(fontsize=8)
        w=windows[key];fig.suptitle(f"Loa 0.20 ≤ z < {w['zhi']:.2f} · log₁₀ M★ ≥ {w['logmass']:.2f}\nRed: 0.75 < g−r < 1.5; blue: g−r < 0.75 · inferred environment probabilities\nSolid: mass/z control within each colour (mass panel: z only); dotted: raw · 95% sky-block intervals",fontsize=12)
        fig.text(.5,.01,'Red and blue controlled curves have different retained support. No robust common red/blue + hard/soft comparison passed the predeclared gate.',ha='center',fontsize=9)
        fig.tight_layout(rect=[0,.035,1,.89]);fig.savefig(OUT/f'properties_{key}.png',dpi=180);fig.savefig(OUT/f'properties_{key}.pdf');plt.close(fig)

def main():
    import fitsio
    from numpy.lib.recfunctions import append_fields
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):raise RuntimeError('Compute allocation required')
    OUT.mkdir(exist_ok=False);allwindows=json.loads((AUDIT/'window_counts.json').read_text())['windows'];windows={k:allwindows[k] for k in ['z020030','z020022']}
    a=np.load(SOURCE/'analysis_input.npy',allow_pickle=False);m=np.log10(a['MASS_CG_15']);results={}
    for key,w in windows.items():
        b=a[(a['Z']>=.2)&(a['Z']<w['zhi'])&(m>=w['logmass'])];assert len(b)==w['rows']
        results[key]=analyse(b)
        if key=='z020022':
            cm=masks(b['GR_OBS']);flag=np.full(len(b),-1,dtype='i2');flag[cm['blue']]=0;flag[cm['red']]=1
            sub=append_fields(b,'COLOUR_CLASS',flag,usemask=False)
            sub=append_fields(sub,'COLOUR_LOGM_LIMIT',np.full(len(b),w['logmass']),usemask=False)
            fitsio.write(OUT/f'LOA_COLOUR_MASS_COMPLETE_{key}.fits',sub,header={'MCFIX':w['logmass'],'MCZLO':.2,'MCZHI':w['zhi'],'COLRED':'0.75 < g-r < 1.5','COLBLUE':'g-r < 0.75'})
            ix=np.linspace(0,len(b)-1,257,dtype=int);fresh=fitsio.read(OUT/f'LOA_COLOUR_MASS_COMPLETE_{key}.fits',rows=ix)
            for col in b.dtype.names:
                eq=np.array_equal(fresh[col],b[col][ix],equal_nan=True) if fresh[col].dtype.kind in 'fc' else np.array_equal(fresh[col],b[col][ix]);assert eq,col
        print('Completed',key,flush=True)
    write_json(OUT/'relationships.json',results);write_json(OUT/'windows.json',windows);render(results,windows)
    write_json(OUT/'MANIFEST.json',{'job':os.environ['SLURM_JOB_ID'],'host':socket.gethostname(),'source_manifest_sha256':sha(SOURCE/'MANIFEST.json'),'colour_audit_manifest_sha256':sha(AUDIT/'MANIFEST.json'),'code':{str(p):sha(p) for p in [Path(__file__),Path(__file__).with_name('run.py'),Path(__file__).with_name('colour_split.py')]},'outputs':{p.name:sha(p) for p in OUT.iterdir() if p.is_file()},'validation':'Every raw and standardized mean independently reconstructed; all source columns roundtripped at 257 rows','interpretation':'No common-colour/hard-soft support passed. Narrow sample is exploratory; soft-only within-colour control has different supports. Bootstrap conditional on fixed selection and predictions.'})
if __name__=='__main__':main()
