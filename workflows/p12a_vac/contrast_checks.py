"""Frozen-product mock attenuation and matched CIGALE comparison; no fitting of VAC."""
import argparse,json,os
from pathlib import Path
import numpy as np
import pandas as pd
import fitsio,healpy as hp
from scipy.spatial import cKDTree
from scipy.stats import rankdata,spearmanr
from astropy.cosmology import Planck18
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from property_environment import digest,VAC,CIG,NAMES,COLORS
ROOT=Path('/pscratch/sd/d/dkololgi/abacus/p10_multiphase/p12a_halo48_candidate_20260924_v1')

def softplus_eigen(t):
 e=np.empty_like(t);e[...,0]=t[...,0];e[...,1]=e[...,0]+np.logaddexp(0,t[...,1]);e[...,2]=e[...,1]+np.logaddexp(0,t[...,2]);return e

def joined_estimates(y,weights,cells,regions,minimum=20,nboot=128):
 """All named estimators share rows, strata, reference mix and bootstrap draws."""
 if not len(y):return None
 _,cells=np.unique(cells,return_inverse=True);_,regions=np.unique(regions,return_inverse=True)
 nc=cells.max()+1;nr=regions.max()+1;flat=cells*nr+regions;counts=[];sums=[];common=np.ones(nc,bool)
 for p in weights.values():
  cnt=np.stack([np.bincount(flat,weights=p[:,j],minlength=nc*nr).reshape(nc,nr) for j in range(4)])
  sy=np.stack([np.bincount(flat,weights=p[:,j]*y,minlength=nc*nr).reshape(nc,nr) for j in range(4)])
  sq=np.stack([np.bincount(cells,weights=p[:,j]**2,minlength=nc) for j in range(4)])
  n=cnt.sum(axis=2);common &= np.all(n*n/np.maximum(sq,1e-30)>=minimum,axis=0)
  counts.append(cnt);sums.append(sy)
 if not common.any():return None
 ref=np.bincount(cells,minlength=nc)[common].astype(float);ref/=ref.sum();C=np.stack(counts)[:,:,common,:];S=np.stack(sums)[:,:,common,:]
 def calc(cc,ss):return (ss/cc*ref).sum(axis=-1)
 means=calc(C.sum(axis=-1),S.sum(axis=-1));rng=np.random.default_rng(28926);bs=[]
 for _ in range(nboot):
  mult=np.bincount(rng.integers(0,nr,nr),minlength=nr);c=C@mult;s=S@mult
  if np.any(c<=0):continue
  bs.append(calc(c,s))
 if len(bs)<nboot*.8:return None
 bs=np.array(bs);dc=bs[:,:,3]-bs[:,:,0]
 return {'rows_input':len(y),'rows_common':int(common[cells].sum()),'cells_common':int(common.sum()),'sky_blocks':int(nr),'bootstrap_draws':len(bs),'estimates':{name:{'means':means[i].tolist(),'contrast':float(means[i,3]-means[i,0]),'interval':np.quantile(dc[:,i],[.16,.84]).tolist()} for i,name in enumerate(weights)},'paired_contrast_differences':{f'{b} minus {a}':{'value':float((means[j,3]-means[j,0])-(means[i,3]-means[i,0])),'interval':np.quantile(dc[:,j]-dc[:,i],[.16,.84]).tolist()} for i,a in enumerate(weights) for j,b in enumerate(weights) if j>i}}

def mock(out):
 import torch
 ds=ROOT/'dataset';po=ROOT/'posterior';au=po/'calibration_audit';marker=ds/'P12A_DATASET_READY.json';ready=json.loads(marker.read_text());complete=json.loads((po/'P12A_COMPLETE.json').read_text())
 assert not ready.get('sealed_phase_opened') and not complete.get('sealed_phase_opened');assert digest(marker)==complete['dataset_marker_sha256'];assert Path(ready['validation']['path']).name=='ph006_selection_sample.npz'
 checkpoint=torch.load(po/'fmpe_estimator.pt',map_location='cpu',weights_only=False);assert checkpoint['dataset_marker_sha256']==digest(marker)
 f=np.load(ds/'ph006_selection_sample.npz');ix=np.load(au/'evaluation_index.npy');draw=np.load(au/'evaluation_samples_scaled.npy',mmap_mode='r');assert draw.shape==(len(ix),512,3)
 truth=f['truth_eigenvalues'][ix];base=f['base_prediction_eigenvalues'][ix];ctx=f['context'][ix];w=f['natural_weight'][ix].astype(float);z=ctx[:,3]
 assert np.allclose(softplus_eigen(f['theta_softplus'][ix]),truth,atol=2e-5)
 p=np.zeros((len(ix),4));coverage=np.zeros((len(ix),3),bool)
 for lo in range(0,len(ix),1000):
  a=softplus_eigen(np.array(draw[lo:lo+1000],dtype='f8')*np.asarray(checkpoint['theta_std'])+np.asarray(checkpoint['theta_mean']));classes=(a>.2).sum(axis=2)
  p[lo:lo+len(a)]=np.stack([(classes==j).mean(axis=1) for j in range(4)],axis=1)
  bounds=np.quantile(a,[.16,.84],axis=1);coverage[lo:lo+len(a)]=(truth[lo:lo+len(a)]>=bounds[0])&(truth[lo:lo+len(a)]<=bounds[1])
 t=(truth>.2).sum(axis=1);oh=np.eye(4)[t];hard=np.eye(4)[p.argmax(axis=1)]
 blockraw=np.column_stack([f['cap'][ix],f['superblock_id'][ix]]);_,blocks=np.unique(blockraw,axis=0,return_inverse=True)
 def summary(k):
  pp=p[k];tt=t[k];ww=w[k];yy=oh[k];hh=hard[k];bb=blocks[k];fold=bb%2
  prev=np.average(yy,axis=0,weights=ww);baseline=np.zeros_like(pp)
  for j in [0,1]:baseline[fold==j]=np.average(yy[fold!=j],axis=0,weights=ww[fold!=j])
  brier=np.average((pp-yy)**2,axis=0,weights=ww);brier0=np.average((baseline-yy)**2,axis=0,weights=ww)
  confusion=np.stack([np.average(hh[tt==j],axis=0,weights=ww[tt==j]) for j in range(4)])
  reli=[]
  for j in range(4):
   bins=[]
   for low,high in zip(np.linspace(0,1,11)[:-1],np.linspace(0,1,11)[1:]):
    keep=(pp[:,j]>=low)&(pp[:,j]<(high if high<1 else 1.0001))
    if keep.sum():bins.append({'predicted':float(np.average(pp[keep,j],weights=ww[keep])),'observed':float(np.average(yy[keep,j],weights=ww[keep])),'rows':int(keep.sum())})
   reli.append(bins)
  # Synthetic unit relation y=true_class/3; not SFR or an HOD physical prediction.
  prop=tt/3;transfer={};rng=np.random.default_rng(672)
  for name,assignment in [('truth',yy),('soft',pp),('argmax',hh)]:
   means=np.sum(assignment*ww[:,None]*prop[:,None],axis=0)/np.sum(assignment*ww[:,None],axis=0)
   transfer[name]={'means':means.tolist(),'contrast':float(means[3]-means[0])}
  null=rng.normal(size=len(tt));nullmeans=np.sum(pp*ww[:,None]*null[:,None],axis=0)/np.sum(pp*ww[:,None],axis=0)
  # Spatial paired resampling of contrast retention; preserve natural weights.
  _,bi=np.unique(bb,return_inverse=True);nb=bi.max()+1;ret=[]
  for _ in range(128):
   m=np.bincount(rng.integers(0,nb,nb),minlength=nb);bw=ww*m[bi];vals=[]
   for assignment in [pp,hh]:
    den=(assignment*bw[:,None]).sum(axis=0);means=(assignment*bw[:,None]*prop[:,None]).sum(axis=0)/den;vals.append(means[3]-means[0])
   ret.append(vals)
  for j,name in enumerate(['soft','argmax']):transfer[name]['interval']=np.quantile(np.array(ret)[:,j],[.16,.84]).tolist()
  return {'rows':int(k.sum()),'class_prevalence':prev.tolist(),'brier_per_class':brier.tolist(),'baseline_brier_per_class':brier0.tolist(),'brier_skill_per_class':(1-brier/brier0).tolist(),'multiclass_brier':float(brier.sum()),'baseline_multiclass_brier':float(brier0.sum()),'logscore_jeffreys512':float(np.average(-np.log((pp[np.arange(len(tt)),tt]*512+.5)/514),weights=ww)),'baseline_logscore':float(np.average(-np.log(np.maximum(baseline[np.arange(len(tt)),tt],1e-12)),weights=ww)),'accuracy':float(np.average(pp.argmax(axis=1)==tt,weights=ww)),'base_eigenvalue_accuracy':float(np.average((base[k]>.2).sum(axis=1)==tt,weights=ww)),'confusion_truth_rows_predicted_columns':confusion.tolist(),'reliability':reli,'synthetic_unit_contrast':transfer,'synthetic_independent_null_soft_contrast':float(nullmeans[3]-nullmeans[0]),'spatial_blocks':int(nb)}
 report={'scope':'Already-exposed ph006 50k evaluation only; natural-weighted, no fresh phase access','synthetic':'Expected proxy y=true_class/3, unit knot-void truth contrast. No physical SFR claim. Independent Gaussian null seed672.','baseline':'Prevalence predicted using opposite parity spatial blocks; no environment model retrained','all':summary(np.ones(len(t),bool)),'wedge_redshift':summary((z>=.2)&(z<.3)),'coverage68':np.average(coverage,axis=0,weights=w).tolist(),'source_hashes':{str(path):digest(path) for path in [marker,ds/'ph006_selection_sample.npz',po/'fmpe_estimator.pt',au/'evaluation_index.npy',au/'evaluation_samples_scaled.npy']}}
 (out/'MOCK.json').write_text(json.dumps(report,indent=2)+'\n');print('MOCK',report['wedge_redshift']['synthetic_unit_contrast'],flush=True)
 return report

def observed(out,cache):
 vh=digest(VAC);assert vh=='cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc'
 cols=['TARGETID','ra','dec','z','MASS_CG','SFR_CG','P_void','P_wall','P_filament','P_cluster'];c=pd.read_parquet(CIG,columns=cols);assert c.TARGETID.is_unique
 ids=fitsio.read(VAC,columns=['TARGETID'])['TARGETID'];order=np.argsort(ids);sortedids=ids[order];at=np.searchsorted(sortedids,c.TARGETID.to_numpy());match=(at<len(ids))&(sortedids[np.minimum(at,len(ids)-1)]==c.TARGETID.to_numpy());cr=np.flatnonzero(match);vr=order[at[cr]];s=np.argsort(vr);c=c.iloc[cr[s]].reset_index(drop=True);vr=vr[s]
 d=fitsio.read(VAC,rows=vr,columns=['TARGETID','RA','DEC','Z','SUPPORTED','P_WEB','QUALITY','BOUNDARY_MPC','NTILDE_MPC3']);assert np.array_equal(d['TARGETID'],c.TARGETID)
 dz=np.abs(c.z.to_numpy()-d['Z'])/(1+d['Z']);sep=np.hypot(((c.ra.to_numpy()-d['RA']+180)%360-180)*np.cos(np.deg2rad(d['DEC'])),c.dec.to_numpy()-d['DEC'])*3600
 mass=c.MASS_CG.to_numpy();sf=c.SFR_CG.to_numpy();ok=d['SUPPORTED']&(dz<.001)&(sep<1)&(d['RA']>=120)&(d['RA']<160)&(d['DEC']>=14.5)&(d['DEC']<30.6)&(d['Z']>=.2)&(d['Z']<.3)&np.isfinite(mass)&np.isfinite(sf)&(mass>0)&(sf>0)
 c=c.loc[ok].reset_index(drop=True);d=d[ok];lm=np.log10(c.MASS_CG.to_numpy());ss=np.log10(c.SFR_CG.to_numpy())-lm;low=(ss<-11).astype(float);eligible=(lm>=8)&(lm<13)
 old=c[['P_void','P_wall','P_filament','P_cluster']].to_numpy(dtype=float);cur=d['P_WEB'].astype(float);assert np.allclose(old.sum(axis=1),1,atol=1e-5) and np.isfinite(old).all();assert np.all(old>=0)
 methods={'current_soft':cur,'current_argmax':np.eye(4)[cur.argmax(axis=1)],'old_soft':old,'old_argmax':np.eye(4)[old.argmax(axis=1)]}
 region=hp.ang2pix(8,d['RA'],d['DEC'],lonlat=True);z=d['Z']
 def compare(mask,dzbin,dmbin,weights=methods):
  k=eligible&mask;nc=int(np.ceil(.1/dzbin))+1;cells=np.floor((z[k]-.2)/dzbin).astype(int)+nc*np.floor((lm[k]-8)/dmbin).astype(int)
  return {name:joined_estimates(y[k],{n:v[k] for n,v in weights.items()},cells,region[k]) for name,y in [('logssfr',ss),('low_ssfr',low)]}
 standard=compare(np.ones(len(d),bool),.025,.25);fine=compare(np.ones(len(d),bool),.01,.1);confident=compare(cur.max(axis=1)>=.8,.01,.1)
 # Direct tracer counts in a 7 Mpc/h sphere, using surrounding full-VAC tracers.
 radius=7/.6766;zt=np.linspace(.15,.55,10001);rt=Planck18.comoving_distance(zt).value
 def xyz(ra,dec,z):
  r=np.interp(z,zt,rt);aa=np.deg2rad(ra);dd=np.deg2rad(dec);return np.column_stack([r*np.cos(dd)*np.cos(aa),r*np.cos(dd)*np.sin(aa),r*np.sin(dd)])
 x=xyz(d['RA'],d['DEC'],z);lo=x.min(axis=0)-radius-.001;hi=x.max(axis=0)+radius+.001;tracers=[]
 with fitsio.FITS(VAC) as f:
  n=f[1].get_nrows()
  for start in range(0,n,200000):
   a=f[1].read(rows=np.arange(start,min(n,start+200000)),columns=['RA','DEC','Z']);xx=xyz(a['RA'],a['DEC'],a['Z']);inside=np.all((xx>=lo)&(xx<=hi),axis=1);tracers.append(xx[inside])
 points=np.concatenate(tracers);tree=cKDTree(points);counts=tree.query_ball_point(x,radius,return_length=True,workers=8)-1;assert np.all(counts>=0)
 expected=4*np.pi*radius**3/3*d['NTILDE_MPC3'];rho=counts/expected
 interior=(d['BOUNDARY_MPC']>radius)&np.isfinite(rho)&(expected>0)&eligible
 ranks=np.full(len(d),-1);density_edges=[]
 # Four density rank groups per narrow z bin, never called web classes; ties stay together.
 zbin=np.floor((z-.2)/.01).astype(int)
 for b in np.unique(zbin[interior]):
  k=interior&(zbin==b);edges=np.quantile(rho[k],[.25,.5,.75]);ranks[k]=np.searchsorted(edges,rho[k],side='right');density_edges.append({'zbin':int(b),'edges':edges.tolist()})
 dp=np.zeros_like(cur);dp[np.flatnonzero(interior),ranks[interior]]=1
 refweights={'current_soft':cur,'current_argmax':methods['current_argmax'],'density_quartiles':dp};reference=compare(interior,.01,.1,refweights)
 # Conditional permutation check: permute current probability vectors within mass/z cells, retaining dependencies within each vector.
 k=eligible;cells=np.floor((z[k]-.2)/.01).astype(int)+11*np.floor((lm[k]-8)/.1).astype(int);_,ci=np.unique(cells,return_inverse=True);rng=np.random.default_rng(717);perm=np.arange(k.sum())
 for cell in np.unique(ci):
  ii=np.flatnonzero(ci==cell);perm[ii]=rng.permutation(ii)
 null=joined_estimates(low[k],{'actual':cur[k],'within_cell_permutation':cur[k][perm]},cells,region[k])
 cache.mkdir(parents=True,exist_ok=False);np.savez(cache/'wedge_checks.npz',TARGETID=d['TARGETID'],mass=lm,logssfr=ss,z=z,current=cur,old=old,region=region,density=rho,counts=counts,interior=interior,density_group=ranks)
 report={'rows':len(d),'vac_sha256':vh,'cache_sha256':digest(CIG),'standard':standard,'fine':fine,'same_current_confident_rows':confident,'density_reference':reference,'permutation':null,'density':{'radius_Mpc':radius,'tracers_in_padded_bbox':len(points),'query_rows':len(d),'interior_rows':int(interior.sum()),'method':'(neighbour count excluding self)/(4pi R^3/3 * frozen NTILDE_MPC3); redshift-space top-hat, not Gaussian tidal environment. No new selection or model fit. Four density ranks within dz=.01; ties not split. Same data, estimator independent of neural weights, not independent observations.','rank_edges':density_edges,'spearman_expected_class_vs_density':float(spearmanr(cur[interior]@np.arange(4),rho[interior]).statistic)},'joined_artifact':str(cache/'wedge_checks.npz'),'joined_sha256':digest(cache/'wedge_checks.npz'),'old_model_provenance':'June 2026 desi_wedge_flowjax_linear_si; cached probabilities averaged over duplicate graph nodes by historical join; lambda threshold .2. Not assumed representative of every earlier model.','vac_unchanged':digest(VAC)==vh}
 assert report['vac_unchanged'];(out/'OBSERVATIONS.json').write_text(json.dumps(report,indent=2)+'\n');print('OBSERVATIONS',fine,flush=True);return report

def figures(out,m,o):
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15});pdf=PdfPages(out/'contrast_checks.pdf')
 def save(fig,name):
  fig.text(.5,.015,'Frozen existing products · mock diagnostics are not DESI truth · observational associations are provisional',ha='center',fontsize=9);fig.savefig(out/(name+'.png'),dpi=170);pdf.savefig(fig);plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(16,5.7));fig.subplots_adjust(top=.79,bottom=.23,wspace=.35);fig.suptitle('Existing mock: discrimination and contrast retention',fontsize=20)
 r=m['wedge_redshift'];im=axs[0].imshow(r['confusion_truth_rows_predicted_columns'],vmin=0,vmax=1,cmap='Blues');axs[0].set_xticks(range(4),NAMES,rotation=25);axs[0].set_yticks(range(4),NAMES);axs[0].set(xlabel='Predicted argmax',ylabel='True class',title='ph006; 0.20 ≤ z < 0.30')
 for i in range(4):
  for j in range(4):axs[0].text(j,i,f"{r['confusion_truth_rows_predicted_columns'][i][j]:.2f}",ha='center',color='white' if r['confusion_truth_rows_predicted_columns'][i][j]>.5 else 'black')
 for j in range(4):axs[1].plot([v['predicted'] for v in r['reliability'][j]],[v['observed'] for v in r['reliability'][j]],marker=['o','s','^','D'][j],label=NAMES[j],color=COLORS[j])
 axs[1].plot([0,1],[0,1],':',color='grey');axs[1].set(xlabel='Predicted probability',ylabel='Observed class frequency',title='Reliability');axs[1].legend(fontsize=9)
 for x,key in enumerate(['truth','soft','argmax']):
  rr=r['synthetic_unit_contrast'][key];axs[2].plot(x,rr['contrast'],'o',color='#2585bb');
  if 'interval' in rr:axs[2].vlines(x,*rr['interval'],color='#2585bb')
 axs[2].set_xticks(range(3),['Truth','Soft weights','Argmax']);axs[2].set(ylim=(0,1.1),ylabel='Recovered / unit true contrast',title='Synthetic proxy; not physical SFR');save(fig,'01_mock_discrimination')
 fig,axs=plt.subplots(1,2,figsize=(13,6));fig.subplots_adjust(left=.2,right=.97,top=.8,bottom=.22,wspace=.35);fig.suptitle('Same galaxies, same fine mass/redshift controls',fontsize=20)
 labels=['current_soft','current_argmax','old_soft','old_argmax']
 for a,metric in zip(axs,['logssfr','low_ssfr']):
  r=o['fine'][metric]
  if r:
   for i,n in enumerate(labels):
    rr=r['estimates'][n];a.hlines(i,*rr['interval'],color='#2585bb');a.plot(rr['contrast'],i,'o',color='#2585bb')
   a.set_title(f"{r['rows_common']:,} shared rows; {r['cells_common']} strata",fontsize=11)
  a.set_yticks(range(4),['Current probabilities','Current argmax','June probabilities','June argmax']);a.invert_yaxis();a.axvline(0,color='grey',ls=':');a.set_xlabel('Knot − void: '+('mean log sSFR [dex]' if metric=='logssfr' else 'low-sSFR fraction'))
 axs[1].set_yticklabels([]);save(fig,'02_same_sample_models')
 fig,axs=plt.subplots(1,2,figsize=(13,6));fig.subplots_adjust(left=.08,right=.98,top=.76,bottom=.23,wspace=.28);fig.suptitle('Current environments and a direct density reference',fontsize=20)
 fig.text(.5,.865,'Identical interior galaxies and fine control strata · density ranks are not void/sheet/filament/knot labels\nDirect counts use a 7 Mpc/h top-hat sphere, surrounding survey tracers and frozen expected density.',ha='center',fontsize=10)
 for a,metric in zip(axs,['logssfr','low_ssfr']):
  r=o['density_reference'][metric]
  if r:
   for key,col,marker in [('current_soft','#2585bb','o'),('current_argmax','#8954af','^'),('density_quartiles','#aa8310','s')]:a.plot(range(4),r['estimates'][key]['means'],label=key.replace('_',' '),color=col,marker=marker)
   a.set_title(f"{r['rows_common']:,} common rows",fontsize=11)
  a.set_xticks(range(4),['Void / Q1','Sheet / Q2','Filament / Q3','Knot / Q4'],rotation=15);a.set_ylabel('Mean log sSFR [dex]' if metric=='logssfr' else 'Low-sSFR fraction');a.legend(fontsize=9)
 save(fig,'03_density_reference');pdf.close()

def test():
 cells=np.repeat([0,1],400);region=np.arange(800)%13;p=np.vstack([np.tile([.4,.3,.2,.1],(400,1)),np.tile([.1,.2,.3,.4],(400,1))]);y=cells*3
 r=joined_estimates(y,{'a':p,'b':p.copy()},cells,region,minimum=10,nboot=16);assert np.allclose(r['estimates']['a']['means'],1.5);assert r['paired_contrast_differences']['b minus a']['value']==0
 hard=np.tile(np.eye(4),(200,1));y=np.tile(np.arange(4),200)+cells*3;r=joined_estimates(y,{'hard':hard},cells,region,minimum=10,nboot=16);assert np.isclose(r['estimates']['hard']['contrast'],3)
 assert np.allclose(softplus_eigen(np.zeros((1,3))),[[0,np.log(2),2*np.log(2)]])
 print('Common-reference composition, paired identity, true contrast and transform tests pass')

if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--test',action='store_true');a.add_argument('--out',type=Path);a.add_argument('--cache',type=Path);args=a.parse_args()
 if args.test:test()
 else:
  args.out.mkdir(parents=True,exist_ok=False);m=mock(args.out);o=observed(args.out,args.cache);figures(args.out,m,o)
  (args.out/'RUN.json').write_text(json.dumps({'script_sha256':digest(__file__),'job_id':os.environ.get('SLURM_JOB_ID'),'scope':'No new model inference, retraining or phase access; saved ph006 and observational VAC only'},indent=2)+'\n')
  (args.out/'ARTIFACT_SHA256.json').write_text(json.dumps({f.name:digest(f) for f in args.out.iterdir() if f.is_file()},indent=2)+'\n')
