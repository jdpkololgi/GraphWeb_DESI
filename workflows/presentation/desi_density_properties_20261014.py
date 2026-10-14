#!/usr/bin/env python3
"""Current CIGALE properties versus verified cached non-neural point counts.
Single-thread login execution explicitly requested by user; no new spatial index.
"""
from pathlib import Path
import argparse, hashlib, json
import numpy as np
import fitsio
from scipy.stats import rankdata, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
ROOT=Path('/global/homes/d/dkololgi/GraphWeb_DESI')
OBS=ROOT/'docs/figures/p12a_contrast_checks_20260928/OBSERVATIONS.json'
PRODUCT=ROOT/'docs/figures/p12a_loa_cigale_20260929/PRODUCT.json'
METHODS=['density','raw_counts','trace','soft_web']
LABELS=['Point density','Raw counts','Inferred trace','Soft web classes']
COLORS=['#F2C14E','#82D9CE','#C28CF5','#E85D9E']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def quartiles(x,zb):
 g=np.zeros(len(x),dtype=int);edges={}
 for b in np.unique(zb):
  k=zb==b;e=np.quantile(x[k],[.25,.5,.75]);g[k]=np.searchsorted(e,x[k],side='right');edges[str(b)]=e.tolist()
 return np.eye(4)[g],edges
def estimate(y,weights,cells,region,x,minimum=20,nboot=128):
 """Identical rows, common cells, reference population and paired bootstrap."""
 ids,c=np.unique(cells,return_inverse=True);_,r=np.unique(region,return_inverse=True)
 nc=len(ids);nr=r.max()+1;flat=c*nr+r;C=[];S=[];X=[];common=np.ones(nc,bool)
 for w in weights.values():
  a=np.stack([np.bincount(flat,weights=w[:,j],minlength=nc*nr).reshape(nc,nr) for j in range(4)])
  b=np.stack([np.bincount(flat,weights=w[:,j]*y,minlength=nc*nr).reshape(nc,nr) for j in range(4)])
  xx=np.stack([np.bincount(flat,weights=w[:,j]*x,minlength=nc*nr).reshape(nc,nr) for j in range(4)])
  q=np.stack([np.bincount(c,weights=w[:,j]**2,minlength=nc) for j in range(4)])
  common &= ((a.sum(2)**2/np.maximum(q,1e-30))>=minimum).all(0)
  C.append(a);S.append(b);X.append(xx)
 if not common.any():raise ValueError('No shared support')
 ref=np.bincount(c,minlength=nc)[common].astype(float);ref/=ref.sum()
 C=np.stack(C)[:,:,common,:];S=np.stack(S)[:,:,common,:];X=np.stack(X)[:,:,common,:]
 def calc(a,b):return (b/a*ref).sum(-1)
 mean=calc(C.sum(-1),S.sum(-1));xm=calc(C.sum(-1),X.sum(-1));rng=np.random.default_rng(20261014);bs=[]
 for _ in range(nboot):
  mult=np.bincount(rng.integers(nr,size=nr),minlength=nr);cc=C@mult
  if (cc>0).all():bs.append(calc(cc,S@mult))
 assert len(bs)>=.8*nboot
 bs=np.array(bs);contr=bs[:,:,-1]-bs[:,:,0]
 out={'input_rows':len(y),'rows':int(common[c].sum()),'cells':int(common.sum()),'sky_blocks':int(nr),'bootstrap_valid':len(bs),'common_cell_ids':ids[common].tolist(),'reference_weights':ref.tolist(),'methods':{}}
 for j,name in enumerate(weights):out['methods'][name]={'mean':mean[j].tolist(),'interval':np.quantile(bs[:,j],[.16,.84],axis=0).tolist(),'mean_log10_density':xm[j].tolist(),'contrast':float(mean[j,-1]-mean[j,0]),'contrast_interval':np.quantile(contr[:,j],[.16,.84]).tolist()}
 out['paired_density_minus_soft_contrast']={'value':float((mean[0,-1]-mean[0,0])-(mean[-1,-1]-mean[-1,0])),'interval':np.quantile(contr[:,0]-contr[:,-1],[.16,.84]).tolist()}
 return out,common[c]
def conditional_rank_corr(x,y,c):
 """Centered within-cell rank correlation; not a partial-Spearman p-value."""
 rx=np.zeros(len(x));ry=rx.copy()
 for b in np.unique(c):
  k=c==b;rx[k]=rankdata(x[k])/k.sum();ry[k]=rankdata(y[k])/k.sum();rx[k]-=rx[k].mean();ry[k]-=ry[k].mean()
 return float(np.corrcoef(rx,ry)[0,1])
def prepare(out):
 o=json.loads(OBS.read_text());cache=Path(o['joined_artifact']);assert sha(cache)==o['joined_sha256']
 a=np.load(cache);receipt=json.loads(PRODUCT.read_text());src=Path(receipt['path'])
 with fitsio.FITS(src) as f:
  assert f[1].get_nrows()==receipt['rows'];ids=f[1].read(columns=['TARGETID'])['TARGETID']
  order=np.argsort(ids);at=np.searchsorted(ids[order],a['TARGETID']);assert np.all(at<len(ids));rows=order[at]
  assert np.array_equal(ids[rows],a['TARGETID']);sort=np.argsort(rows);undo=np.argsort(sort)
  cols=['TARGETID','Z','SUPPORTED','P_WEB','EIGENVALUE_MEAN','CIGALE_MATCHED','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15','AGE_CG_15','AV_CG_15','FAGN_CG_15']
  d=f[1].read(rows=rows[sort],columns=cols)[undo]
 assert np.array_equal(d['TARGETID'],a['TARGETID']) and np.allclose(d['P_WEB'],a['current'],atol=1e-7) and np.allclose(d['Z'],a['z'],atol=1e-9)
 assert np.all(np.diff(d['EIGENVALUE_MEAN'],axis=1)>=-1e-6)
 mass=d['MASS_CG_15'];sfr=d['SFR_CG_15'];age=d['AGE_CG_15'];av=d['AV_CG_15'];agn=d['FAGN_CG_15']
 valid=a['interior']&d['SUPPORTED']&d['CIGALE_MATCHED']&(np.char.strip(d['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')
 for v in [mass,sfr,age]:valid &= np.isfinite(v)&(v>0)
 valid &= np.isfinite(av)&(av>=0)&np.isfinite(agn)&(agn>=0)&(agn<=1)&np.isfinite(a['density'])&(a['density']>0)
 k=np.flatnonzero(valid);d=d[k];lm=np.log10(mass[k]);ls=np.log10(sfr[k]);ss=ls-lm;rho=a['density'][k];cnt=a['counts'][k];tr=d['EIGENVALUE_MEAN'].sum(1).astype(float);zb=np.floor((d['Z']-.2+1e-10)/.01).astype(int);cells=zb+11*np.floor((lm-8+1e-10)/.1).astype(int);region=a['region'][k]
 weights={};edges={}
 for name,x in [('density',rho),('raw_counts',cnt),('trace',tr)]:weights[name],edges[name]=quartiles(x,zb)
 weights['soft_web']=d['P_WEB'].astype(float)
 properties={'mass':(lm,'log₁₀ stellar mass [M☉]'),'sfr':(ls,'log₁₀ SFR [M☉ yr⁻¹]'),'ssfr':(ss,'log₁₀ sSFR [yr⁻¹]'),'quenched':((ss<-11).astype(float),'Quenched fraction'),'age':(np.log10(age[k]/1000),'log₁₀ fitted age [Gyr]'),'av':(av[k],'Dust attenuation Aᵥ [mag]'),'fagn':(agn[k],'AGN fraction')}
 stats={}
 for name,(y,label) in properties.items():
  c=zb if name=='mass' else cells;r,keep=estimate(y,weights,c,region,np.log10(rho));r['label']=label;r['control']='redshift only' if name=='mass' else 'mass + redshift';r['correlations']={}
  for key,x in [('density',rho),('raw_counts',cnt),('trace',tr)]:r['correlations'][key]={'spearman':float(spearmanr(x[keep],y[keep]).statistic),'within_cell_rank':conditional_rank_corr(x[keep],y[keep],c[keep])}
  stats[name]=r;print(name,'rows',r['rows'],'density',r['methods']['density']['contrast'],'soft',r['methods']['soft_web']['contrast'],flush=True)
 np.savez_compressed(out/'matched_galaxies.npz',TARGETID=d['TARGETID'],catalogue_rows=rows[k],z=d['Z'],mass=lm,log_sfr=ls,log_ssfr=ss,log_age_Gyr=properties['age'][0],av=av[k],fagn=agn[k],density=rho,counts=cnt,trace=tr,P_WEB=d['P_WEB'],region=region,mass_z_cells=cells,z_cells=zb,**{f'weights_{n}':w for n,w in weights.items()})
 dump(out/'statistics.json',stats)
 dump(out/'provenance.json',{'source_receipt':receipt,'source_receipt_sha256':sha(PRODUCT),'full_catalogue_sha256_recomputed':False,'binding':'Catalogue row count, exact TARGETID and cached current P_WEB/redshift parity checked. Full catalogue receipt hash inherited; extracted rows hashed below.','density_cache':str(cache),'density_cache_sha256':sha(cache),'extracted_rows_sha256':sha(out/'matched_galaxies.npz'),'script_sha256':sha(__file__),'cached_queries':len(a['TARGETID']),'common_valid_properties':len(d),'scope':'Historical density-cache wedge RA120-160 Dec14.5-30.6 z.20-.30; current 15-band CIGALE join. Historical cache required positive old mass/SFR; this remains a sample selection, not a full-footprint or volume-complete analysis.','density':o['density'],'property_selection':'All seven properties share finite physically valid input rows; mass, SFR, age >0; AV>=0; FAGN in [0,1]; matched galaxy and supported interior; density>0. Age is the catalogue SED-fit AGE parameter, not verified mass-weighted age.','quartile_edges':edges,'controls':'Quartiles defined within dz=.01 before property analysis; ties unsplit. Shared mass(.1dex)/z(.01) strata across all four estimators, Kish>=20 in every group. Mass property uses z-only control. Row count differs between mass and other properties. All other panels use identical galaxies and reference mix.','uncertainty':'128 paired nside8 sky-block resamples;16-84% intervals; fixed strata/reference/quartile thresholds. Spatial-sampling uncertainty only.','correlation':'Spearman and centered within-cell rank correlation on common rows; latter is descriptive, no p-values. Quenched correlation is with a binary indicator.','execution':'User explicitly requested login-node execution. Single thread, cached counts, no Slurm action.','no_new_inference':True})
def render(out):
 r=json.loads((out/'statistics.json').read_text())
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'figure.facecolor':'#080C13','axes.facecolor':'#080C13','text.color':'#EDF1F5','axes.labelcolor':'#EDF1F5','xtick.color':'#BBC5D2','ytick.color':'#BBC5D2','axes.edgecolor':'#8995A5','pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'#080C13'})
 for mode in ['point_density_properties','environment_comparison']:
  f,axs=plt.subplots(2,4,figsize=(18,10));f.subplots_adjust(left=.065,right=.98,bottom=.15,top=.79,wspace=.36,hspace=.48)
  f.suptitle('CIGALE properties versus independent point density' if mode=='point_density_properties' else 'Independent density and inferred environments: matched comparison',fontsize=23,y=.95)
  f.text(.5,.898,'Current 15-band CIGALE · cached DESI wedge, 0.20 ≤ z < 0.30 · 7 h⁻¹ Mpc neighbour counts',ha='center',fontsize=14)
  for ax,(name,v) in zip(axs.flat,r.items()):
   for j,method in enumerate(METHODS):
    if mode=='point_density_properties' and method!='density':continue
    d=v['methods'][method];x=np.array(d['mean_log10_density']) if mode=='point_density_properties' else np.arange(4);mu=np.array(d['mean']);ci=np.array(d['interval'])
    ax.plot(x,mu,marker=['o','s','D','^'][j],color=COLORS[j],lw=1.8,label=LABELS[j]);ax.vlines(x,ci[0],ci[1],color=COLORS[j],lw=1.3)
   ax.set_title(v['label'],fontsize=12,pad=12);ax.grid(axis='y',alpha=.15);ax.spines[['top','right']].set_visible(False)
   ax.text(.03,.96,f"{v['rows']:,} galaxies\n{v['control']}",transform=ax.transAxes,va='top',fontsize=9,color='#BBC5D2')
   if name=='quenched':ax.set_ylim(0,1)
   else:ax.margins(y=.30)
   if mode=='point_density_properties':ax.set_xlabel('Mean log₁₀(N / expected N)',fontsize=10)
   else:ax.set_xticks(range(4),['1','2','3','4']);ax.set_xlabel('Environment group',fontsize=10)
  axs.flat[-1].axis('off');axs.flat[-1].text(0,.95,'Same galaxies and control mix\nfor every estimator in each panel.\n\nQuartiles: low → high within Δz=0.01.\nSoft classes: void → wall →\nfilament → cluster.\n\nGroups are different selections;\nmatching amplitudes is a\nconsistency check.\n\nBars: 16–84% sky-block bootstrap.',va='top',linespacing=1.5,fontsize=11)
  if mode=='environment_comparison':f.legend(*axs.flat[0].get_legend_handles_labels(),ncol=4,loc='upper center',bbox_to_anchor=(.5,.865),frameon=False,fontsize=12)
  else:f.text(.5,.85,'Each point is a density quartile; x and y share the same control weights.',ha='center',fontsize=13,color='#F2C14E')
  f.text(.5,.067,'Common 0.1-dex mass / 0.01-redshift strata · stellar mass controls redshift only · quenched: log₁₀ sSFR < −11',ha='center',fontsize=12)
  f.text(.5,.028,'Counts independent of neural weights; same survey and selection · AGE is the fitted CIGALE parameter · provisional inferred environments',ha='center',fontsize=11,color='#BBC5D2')
  for ext in ['png','pdf','svg']:f.savefig(out/f'{mode}.{ext}',dpi=180,bbox_inches='tight')
  plt.close(f);im=Image.open(out/f'{mode}.png');im.thumbnail((1600,1000));im.save(out/f'{mode}_preview.png')
def test():
 cells=np.repeat([0,1],400);regions=np.arange(800)%13;w=np.tile(np.eye(4),(200,1));y=np.tile(np.arange(4),200)+cells*3
 r,k=estimate(y,{n:w for n in METHODS},cells,regions,np.zeros(800),minimum=10,nboot=16)
 assert k.all() and np.isclose(r['methods']['density']['contrast'],3) and r['paired_density_minus_soft_contrast']['value']==0
 assert np.isclose(conditional_rank_corr(y,y,cells),1)
 q,_=quartiles(np.repeat([0,1,2,3],200),np.zeros(800));assert np.array_equal(q.argmax(1),np.repeat([0,1,2,3],200))
 print('Shared-reference, paired-identity, within-cell-rank and tie tests passed')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path);p.add_argument('--test',action='store_true');p.add_argument('--render-only',action='store_true');a=p.parse_args()
 if a.test:test()
 else:
  a.out.mkdir(parents=True,exist_ok=True)
  if not a.render_only:prepare(a.out)
  render(a.out)
