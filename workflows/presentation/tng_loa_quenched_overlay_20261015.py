"""Cross-population rank-aligned quenching curves; no observational truth claim."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from loa_tng_abacus_overlay_20261015 import property_data
H=Path('/global/u2/d/dkololgi');R=H/'fmsu_tokyo_20261010';O=R/'tng_loa_quenched_overlay_20261015';O.mkdir(exist_ok=True)
def quartile(v):
 _,inv,n=np.unique(v,return_inverse=True,return_counts=True);e=np.cumsum(n);s=e-n;N=len(v)
 return np.array([np.maximum(0,np.minimum(e,(k+1)*N/4)-np.maximum(s,k*N/4))/n for k in range(4)]).T[inv]
rng=np.random.default_rng(331);checks=0
# Per mass-bin curves. Loa groups additionally formed within dz=.01.
def strata(y,trace,mb,zb,blocks):
 global checks
 _,blocks=np.unique(blocks,return_inverse=True);nb=blocks.max()+1;draw=rng.multinomial(nb,np.ones(nb)/nb,size=200);out={}
 for b in np.unique(mb):
  point=[];boots=[];sizes=[]
  for z in np.unique(zb[mb==b]):
   ii=np.flatnonzero((mb==b)&(zb==z))
   if len(ii)<80:continue
   w=quartile(trace[ii]);den=np.array([np.bincount(blocks[ii],weights=w[:,k],minlength=nb) for k in range(4)]).T;num=np.array([np.bincount(blocks[ii],weights=w[:,k]*y[ii],minlength=nb) for k in range(4)]).T
   v=num.sum(0)/den.sum(0);assert np.allclose(v,(w.T@y[ii])/w.sum(0),atol=1e-12);checks+=1
   point.append(v);boots.append((draw@num)/(draw@den));sizes.append(len(ii))
  if sizes:
   weights=np.array(sizes)/sum(sizes);out[int(b)]={'n':sum(sizes),'mean':weights@np.array(point),'boot':np.einsum('i,ijk->jk',weights,np.array(boots))}
 return out
T=R/'tng_density_control_20261015';d=np.load(T/'catalogue.npz');t=np.load(T/'truth.npz');f=np.load(T/'features.npz');p=np.load(T/'predictions.npz');ids=f['test_ids'];m=np.log10(d['mass_msun'][ids]);ids=ids[(m>=9)&(m<11.5)];m=np.log10(d['mass_msun'][ids]);mb=np.floor((m-9)/.25).astype(int);q=(d['sfr'][ids]/d['mass_msun'][ids]<1e-11).astype(float);block=np.ravel_multi_index(np.minimum((d['pos_hmpc'][ids]/205*5).astype(int),4).T,(5,5,5));tg={}
for rad in [2,7]:
 for arm in ['truth','dense','thin']:
  names=[None] if arm=='truth' else [f'R{rad}_dense_seed{s}' if arm=='dense' else f'R{rad}_thin_{s}_seed{s}' for s in [20261015,20261016,20261017]]
  # Reset block draws to align spatial resampling across fit seeds and model arms.
  ss=[]
  for name in names:
   rng=np.random.default_rng(331);e=t[f'eig_R{rad}_n512'][ids] if name is None else p[name][ids];ss.append(strata(q,e.sum(1),mb,np.zeros(len(mb),int),block))
  tg[(rad,arm)]={b:{'n':ss[0][b]['n'],'mean':np.mean([s[b]['mean'] for s in ss],0),'boot':np.mean([s[b]['boot'] for s in ss],0)} for b in ss[0]}
a=np.load(R/'loa_tng_abacus_overlay_20261015/loa_predictions.npz');phot=np.load(R/'legacy_selection_tests_20261014/p12a_targeting.npz');assert np.array_equal(a['TARGETID'],phot['TARGETID']);good,y=property_data(a);good&=(phot['r']<phot['rlimit'])&(y['mass']<11.5)
lmb=np.floor((y['mass'][good]-9)/.25+1e-10).astype(int);lzb=np.floor((a['z'][good]-.2)/.01+1e-8).astype(int)
rng=np.random.default_rng(991);lo=strata(y['quenched'][good],a['abacus'][good].sum(1),lmb,lzb,a['region'][good]);results={}
def aggregate(stats,bins,weights):
 v=sum(w*stats[b]['mean'] for b,w in zip(bins,weights));bt=sum(w*stats[b]['boot'] for b,w in zip(bins,weights));return {'mean':v.tolist(),'ci95':np.quantile(bt,[.025,.975],axis=0).tolist(),'contrast_pp':float((v[-1]-v[0])*100),'contrast_ci95_pp':(100*np.quantile(bt[:,-1]-bt[:,0],[.025,.975])).tolist(),'rows':sum(stats[b]['n'] for b in bins),'mass_bins':list(bins),'mass_weights':list(weights)}
for rad in [2,7]:
 fig,axs=plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
 for row,high in enumerate([False,True]):
  for col,matched in enumerate([False,True]):
   ax=axs[row,col];models={'TNG truth':tg[(rad,'truth')],'TNG dense prediction':tg[(rad,'dense')],'TNG sparse, retrained':tg[(rad,'thin')],'P12-A on Loa (R7)':lo};records={}
   shared=sorted(set(lo)&set(tg[(rad,'truth')]));shared=[b for b in shared if not high or b>=6]
   for label,stats in models.items():
    bins=shared if matched else sorted(b for b in stats if not high or b>=6)
    counts=np.array([lo[b]['n'] if matched else stats[b]['n'] for b in bins]);weights=counts/counts.sum();records[label]=aggregate(stats,bins,weights)
   key=f'R{rad}_{"high" if high else "full"}_{"matched" if matched else "native"}';results[key]=records
   styles=[('#222222','o-'),('#1685ad','s--'),('#d76b16','^:'),('#a42289','D-')]
   for (label,rec),(color,style) in zip(records.items(),styles):
    v=np.array(rec['mean'])*100;ci=np.array(rec['ci95'])*100;ax.plot(range(1,5),v,style,color=color,label=label);ax.fill_between(range(1,5),ci[0],ci[1],color=color,alpha=.1)
   ax.set_title(('High mass: 10.5 ≤ log M < 11.5' if high else ('Shared support: 9.5 ≤ log M < 11.5' if matched else 'Full range: 9 ≤ log M < 11.5'))+'\n'+('Same Loa mass-bin mixture' if matched else 'Each sample’s own mass mixture'),fontsize=11)
   ax.set_xticks([1,2,3,4],['Lowest','2','3','Highest']);ax.grid(alpha=.16)
   if col==0:ax.set_ylabel('Stopped forming stars (%)')
   if row==1:ax.set_xlabel('Environment quarter, ranked within mass bins')
 fig.suptitle(f'TNG R{rad} truth/predictions versus P12-A Loa R7: quenched fractions',fontsize=15,y=.995)
 fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.96),ncol=4,fontsize=9,frameon=False)
 fig.text(.5,.012,'Loa: 0.20 ≤ z < 0.30; CG15 properties; ranks also controlled within Δz = 0.01. TNG: z = 0, instantaneous SFR.\nShared threshold sSFR < 10⁻¹¹ yr⁻¹; 0.25-dex mass bins; both samples capped at log M < 11.5.\nShading: 95% spatial intervals. TNG predictions are three-fit geometry-only MLPs, not GAT/SBI.\nEqual rank groups do not imply equal physical environments; matching mass does not match epoch, luminosity or SFR measurement.',ha='center',fontsize=9)
 fig.tight_layout(rect=[0,.14,1,.90]);fig.savefig(O/f'R{rad}_tng_loa_overlay.png',dpi=175,bbox_inches='tight');fig.savefig(O/f'R{rad}_tng_loa_overlay.pdf',bbox_inches='tight');plt.close(fig)
# Check against previous independently aggregated full-population TNG curves.
prior=json.loads((R/'supervisor_environment_figures_20261015/RESULTS.json').read_text())['tng'];parity=[]
for rad in [2,7]:
 for high in ['full','high']:
  for arm,label in [('truth','TNG truth'),('dense','TNG dense prediction'),('thin','TNG sparse, retrained')]:
   old=prior[f'R{rad}_{"all" if high=="full" else "high"}_{arm}']['mean'];new=results[f'R{rad}_{high}_native'][label]['mean'];err=float(np.max(np.abs(np.array(old)[:,0]-new)));assert err<1e-12;parity.append(err)
receipt={'curves':results,'direct_cell_mean_checks':checks,'previous_TNG_curve_checks':len(parity),'maximum_previous_curve_error':max(parity),'loa_eligible_before_cell_support':int(good.sum()),'scope':'fixed cached predictions, no retraining; Loa rank groups within mass/redshift; TNG within mass; shared mass weights equal eligible Loa counts; cells>=80; no claim of equal physical environment or identical SFR measurement','sources':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in [T/'catalogue.npz',T/'truth.npz',T/'predictions.npz',R/'loa_tng_abacus_overlay_20261015/loa_predictions.npz']}}
(O/'RESULTS.json').write_text(json.dumps(receipt,indent=2,allow_nan=False));print('DONE',checks,'checks',flush=True)
for key,records in results.items():print(key,{k:round(v['contrast_pp'],2) for k,v in records.items()},flush=True)
