"""Same-truth TNG mass/density comparison and frozen Loa selection summaries."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path('/global/u2/d/dkololgi');R=H/'fmsu_tokyo_20261010';O=R/'supervisor_environment_figures_20261015';O.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):
 fig.savefig(O/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(O/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def groups(v):
 _,inv,n=np.unique(v,return_inverse=True,return_counts=True);end=np.cumsum(n);start=end-n;total=len(v)
 return np.column_stack([np.maximum(0,np.minimum(end,(k+1)*total/4)-np.maximum(start,k*total/4))/n for k in range(4)])[inv]
d=np.load(R/'tng_density_control_20261015/catalogue.npz');t=np.load(R/'tng_density_control_20261015/truth.npz');f=np.load(R/'tng_density_control_20261015/features.npz');p=np.load(R/'tng_density_control_20261015/predictions.npz')
ids=f['test_ids'];lm=np.log10(d['mass_msun'][ids]);ids=ids[(lm>=9)&(lm<11.5)];lm=np.log10(d['mass_msun'][ids]);bins=np.floor((lm-9)/.25).astype(int)
y=np.column_stack([d['sfr'][ids]/d['mass_msun'][ids]<1e-11,d['g_r'][ids],np.log10(np.maximum(d['sfr'][ids]/d['mass_msun'][ids],1e-14))]).astype(float)
bb=np.ravel_multi_index(np.minimum((d['pos_hmpc'][ids]/205*5).astype(int),4).T,(5,5,5));_,bb=np.unique(bb,return_inverse=True);nb=bb.max()+1;draw=np.random.default_rng(711).multinomial(nb,np.ones(nb)/nb,size=200)
results={};checks=0
for radius in [2,7]:
 for high in [False,True]:
  use=lm>=10.5 if high else np.ones(len(ids),bool)
  for arm in ['truth','dense','thin','frozen_thin']:
   names=[None] if arm=='truth' else [f'R{radius}_dense_seed{s}' if arm=='dense' else f'R{radius}_thin_{s}_seed{s}' if arm=='thin' else f'R{radius}_frozen_thin_seed{s}' for s in [20261015,20261016,20261017]]
   points=[];boots=[]
   for name in names:
    e=t[f'eig_R{radius}_n512'][ids] if name is None else p[name][ids];trace=e.sum(1);point=np.zeros((4,3));boot=np.zeros((200,4,3))
    for cell in np.unique(bins[use]):
     ii=np.flatnonzero(use&(bins==cell));w=groups(trace[ii]);ref=len(ii)/use.sum()
     direct=w.T@y[ii]/w.sum(0)[:,None];counts=np.array([np.bincount(bb[ii],weights=w[:,k],minlength=nb) for k in range(4)]).T
     sums=np.stack([np.array([np.bincount(bb[ii],weights=w[:,k]*y[ii,j],minlength=nb) for k in range(4)]).T for j in range(3)],axis=-1)
     mean=sums.sum(0)/counts.sum(0)[:,None];assert np.allclose(mean,direct,atol=1e-12);checks+=1
     den=draw@counts;bs=np.stack([draw@sums[:,:,j] for j in range(3)],axis=-1)/den[:,:,None]
     point+=ref*mean;boot+=ref*bs
    points.append(point);boots.append(boot)
   point=np.mean(points,0);boot=np.mean(boots,0);results[f'R{radius}_{"high" if high else "all"}_{arm}']={'n':int(use.sum()),'mean':point.tolist(),'ci95':np.quantile(boot,[.025,.975],axis=0).tolist(),'contrast':(point[-1]-point[0]).tolist(),'contrast_ci95':np.quantile(boot[:,-1]-boot[:,0],[.025,.975],axis=0).tolist()}
colors={'truth':'#222222','dense':'#157aab','thin':'#d76b16','frozen_thin':'#9445a0'};labels={'truth':'True environment (no ML)','dense':'Prediction: dense training','thin':'Prediction: sparse, retrained','frozen_thin':'Dense model on sparse inputs'}
for prop,pname in enumerate(['quenched','colour','ssfr']):
 fig,axs=plt.subplots(2,2,figsize=(11,7.8),sharex=True,sharey=True)
 for row,rad in enumerate([2,7]):
  for col,high in enumerate([False,True]):
   ax=axs[row,col]
   for arm in ['truth','dense','thin']:
    rec=results[f'R{rad}_{"high" if high else "all"}_{arm}'];v=np.array(rec['mean'])[:,prop];ci=np.array(rec['ci95'])[:,:,prop];factor=100 if prop==0 else 1
    ax.plot(range(1,5),v*factor,{'truth':'o-','dense':'s--','thin':'^:'}[arm],c=colors[arm],label=labels[arm]);ax.fill_between(range(1,5),ci[0]*factor,ci[1]*factor,color=colors[arm],alpha=.12)
   ax.set_title(f'R = {rad} Mpc/h | '+('Mass > 3.16 × 10¹⁰ M☉' if high else 'Mass > 10⁹ M☉'));ax.set_xticks([1,2,3,4],['Lowest','2','3','Highest']);ax.grid(alpha=.15)
   if col==0:ax.set_ylabel(['Stopped forming stars (%)','Mean g − r colour','Mean log₁₀(sSFR / yr⁻¹)'][prop])
   if row==1:ax.set_xlabel('Environment quarter (within stellar-mass bins)')
 fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.965),ncol=3,frameon=False,fontsize=10)
 fig.suptitle('TNG: changing the galaxy population versus removing neighbours',y=1.01,fontsize=15)
 fig.text(.5,.005,'Same held-out galaxies and matter field. Sparse tracers: 283,969 → 37,545.\nPredictions: geometry-only regressors (not GAT/SBI), mean of 3 fits; shaded 95% spatial intervals.\nMass < 10¹¹⋅⁵ M☉ throughout; 66,691 full / 9,338 high-mass queries. Groups defined separately per model.',ha='center',fontsize=9)
 fig.tight_layout(rect=[0,.09,1,.91]);save(fig,'tng_truth_prediction_'+pname)
# Summary contrast makes mass-selection versus sparsity easiest to compare.
fig,axs=plt.subplots(1,2,figsize=(11,4.7),sharey=True)
for ax,rad in zip(axs,[2,7]):
 for j,arm in enumerate(['truth','dense','thin','frozen_thin']):
  yy=[];low=[];hi=[]
  for high in [False,True]:
   r=results[f'R{rad}_{"high" if high else "all"}_{arm}'];v=r['contrast'][0]*100;yy.append(v);low.append(v-100*r['contrast_ci95'][0][0]);hi.append(100*r['contrast_ci95'][1][0]-v)
  ax.errorbar(np.array([0,1])+(j-1.5)*.045,yy,yerr=[low,hi],fmt='o--' if arm=='frozen_thin' else 'o-',color=colors[arm],label=labels[arm],capsize=3)
 ax.set(title=f'R = {rad} Mpc/h',xticks=[0,1],xticklabels=['Full mass range','High-mass cut']);ax.grid(alpha=.15)
axs[0].set_ylabel('Highest − lowest quarter: quenched difference (points)');fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=2,fontsize=9);fig.tight_layout(rect=[0,.13,1,1]);save(fig,'tng_contrast_summary')
a=json.loads((R/'loa_actual_gat_overlay_20261015/RESULTS.json').read_text());print(a['paired_cuts']['legacy_targetfit_full']['methods'].keys(),flush=True)
fig,axs=plt.subplots(1,3,figsize=(14,4.4));curve=a['curve'];axs[0].plot([(v['zlo']+v['zhi'])/2 for v in curve],[v['logmass90'] for v in curve],'o-',color='#d76b16',label='Estimated 90% completeness limit');axs[0].axhline(9,c='#157aab',ls='--',label='A simple mass > 10⁹ cut');axs[0].set(xlabel='Redshift',ylabel='log₁₀(stellar mass / M☉)',title='A 10⁹ cut is not completeness');axs[0].legend(fontsize=8)
for idx,kind in [(1,'raw'),(2,'paired')]:
 for cut,col,label in [('uncut','#157aab','Mass > 10⁹ parent'),('complete90','#d76b16','90% completeness cut')]:
  if kind=='raw':r=a['scenarios']['legacy_targetfit_full_'+cut]['suites']['soft_raw']['quenched']['methods']['abacus']
  else:r=a['paired_cuts']['legacy_targetfit_full']['methods']['abacus_'+('parent' if cut=='uncut' else 'cut')]
  v=np.array(r['mean'])*100;ci=np.array(r['interval'])*100;axs[idx].plot(range(4),v,'o-',c=col,label=label);axs[idx].fill_between(range(4),ci[0],ci[1],color=col,alpha=.13)
  axs[idx].text(.03,.94-(cut=='complete90')*.1,f'{label}: {r["contrast"]*100:.2f} points',transform=axs[idx].transAxes,color=col,fontsize=9)
 axs[idx].set(xticks=range(4),xticklabels=['Void','Wall','Filament','Cluster'],ylabel='Stopped forming stars (%)',title='Raw population comparison' if kind=='raw' else 'Same mass/redshift mixture');axs[idx].grid(alpha=.15)
fig.suptitle('P12-A Loa, 0.20 ≤ z < 0.30: completeness does not restore the old contrast',fontsize=14,y=1.02)
fig.text(.5,-.015,'Probability-weighted classes; shaded 16–84% spatial intervals. Existing fixed predictions, no retraining.\nRaw parent/cut: 90,952 / 52,620 galaxies. Paired comparison: 67 shared mass/redshift cells; 52,340 parent rows.',ha='center',fontsize=9)
fig.tight_layout();save(fig,'loa_completeness_explanation')
# Full seven-property original data, isolate Abacus for clearer reading.
props=list(a['scenarios']['legacy_targetfit_full_uncut']['suites']['soft_raw'])
fig,axs=plt.subplots(2,4,figsize=(15,7));axs=axs.ravel()
for ax,prop in zip(axs,props):
 for cut,col,label in [('uncut','#157aab','Parent'),('complete90','#d76b16','90% cut')]:
  rr=a['scenarios']['legacy_targetfit_full_'+cut]['suites']['soft_raw'][prop]['methods']['abacus'];v=np.array(rr['mean']);ci=np.array(rr['interval']);factor=100 if prop=='quenched' else 1
  ax.plot(range(4),v*factor,'o-',c=col,label=label);ax.fill_between(range(4),ci[0]*factor,ci[1]*factor,color=col,alpha=.14)
 ax.set(title=prop,xticks=range(4),xticklabels=['Void','Wall','Filament','Cluster']);ax.tick_params(axis='x',labelsize=8);ax.grid(alpha=.15)
axs[-1].axis('off');axs[0].legend();fig.suptitle('P12-A Loa: raw property curves before/after the completeness cut\nPopulation changes between curves; no mass/redshift balancing in these panels');fig.tight_layout(rect=[0,0,1,.91]);save(fig,'loa_all_properties_raw')
(O/'RESULTS.json').write_text(json.dumps({'tng':results,'independent_cell_mean_checks':checks,'method':'Fixed 0.25-dex mass cells; fractional quartile ties; 200 paired spatial draws; mean of three fits; no fit uncertainty. High-mass cut changes queries only, not retraining or tracer population.','loa_source':str(R/'loa_actual_gat_overlay_20261015/RESULTS.json')},indent=2))
for rad,expect in [(2,.494999325246285),(7,.3272855407776163)]:assert abs(results[f'R{rad}_all_truth']['contrast'][0]-expect)<1e-10
(O/'MANIFEST.json').write_text(json.dumps({str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(O.glob('*')) if p.name!='MANIFEST.json'},indent=2))
print('DONE',checks,flush=True)
