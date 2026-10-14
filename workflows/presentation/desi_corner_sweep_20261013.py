#!/usr/bin/env python3
"""Source-bound descriptive eigenvalue mixtures and threshold curves; saved draws only."""
from pathlib import Path
import argparse,json,hashlib,os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter
from PIL import Image
ROOT=Path(os.environ.get('GRAPHWEB_REPO_ROOT','/global/homes/d/dkololgi/GraphWeb_DESI')).resolve()
DECK=Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010')
VACROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1')
NAMES=['Abacus truth','Abacus posterior','DESI posterior'];COLORS=['#C0C0C0','#5B8FF9','#F04D9A']
WEB=['Void','Wall','Filament','Cluster'];WC=['#9CF6D5','#5B8FF9','#E85D9E','#F2C14E'];STYLES=['-','--',':']
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def aggregate(e,w,edges,thresholds):
 assert np.isfinite(e).all() and np.all(np.diff(e,axis=1)>=0)
 den=w.sum();hist=np.array([np.histogram(e[:,j],edges,weights=w)[0]/den/np.diff(edges) for j in range(3)])
 pairs=[(0,1),(0,2),(1,2)];h2=np.array([np.histogram2d(e[:,i],e[:,j],bins=[edges,edges],weights=w)[0]/den for i,j in pairs])
 cdf=[]
 for j in range(3):
  ix=np.searchsorted(thresholds,e[:,j],side='left');c=np.bincount(ix,weights=w,minlength=len(thresholds)+1).cumsum()[:-1]/den;cdf.append(c)
 f1,f2,f3=cdf;frac=np.array([f3,f2-f3,f1-f2,1-f1]).T
 assert np.all(frac>=-1e-10) and np.allclose(frac.sum(1),1)
 outmass=[float(w[(e[:,j]<edges[0])|(e[:,j]>edges[-1])].sum()/den) for j in range(3)]
 return hist,h2,frac,outmass

def prepare(out):
 import torch,fitsio
 torch.set_num_threads(1);sources=[]
 def bind(path,expected=None):
  value=sha(path)
  if expected:assert value==expected,str(path)
  sources.append({'path':str(path),'sha256':value})
 a02=json.loads((DECK/'a02_histograms/provenance_qa.json').read_text());bindings={Path(v['path']).name:v for v in a02['source_bindings']}
 for name in ['evaluation_samples_scaled.npy','evaluation_index.npy','ph006_selection_sample.npz','fmpe_estimator.pt','P12A_DATASET_READY.json']:
  v=bindings[name];bind(v['path'],v['sha256'])
 assert json.loads(Path(bindings['P12A_DATASET_READY.json']['path']).read_text())['validation_phase']=='ph006'
 index=np.load(bindings['evaluation_index.npy']['path'])
 with np.load(bindings['ph006_selection_sample.npz']['path']) as d:
  z=d['context'][index,3];rows=np.flatnonzero((z>=.2)&(z<.3));vi=index[rows];truth=d['truth_eigenvalues'][vi].astype(float);weight=d['natural_weight'][vi].astype(float)
  assert np.isin(d['fold'][vi],[2,3,4]).all() and (d['random_support'][vi]>0).all()
 ck=torch.load(bindings['fmpe_estimator.pt']['path'],map_location='cpu',weights_only=False);raw=np.load(bindings['evaluation_samples_scaled.npy']['path'],mmap_mode='r')
 t=np.asarray(raw[rows]*ck['theta_std']+ck['theta_mean'],dtype=float);e=np.empty_like(t);e[...,0]=t[...,0];e[...,1]=e[...,0]+np.logaddexp(0,t[...,1]);e[...,2]=e[...,1]+np.logaddexp(0,t[...,2]);del t
 m=json.loads((VACROOT/'FULL_VAC_COMPLETE.json').read_text());bind(m['catalogue'],m['catalogue_sha256']);bind(VACROOT/'DRAW_INDEX.json',m['draw_index_sha256']);drawindex=json.loads((VACROOT/'DRAW_INDEX.json').read_text())['shards']
 bind(VACROOT/'SOURCE_MANIFEST.json',m['source_manifest_sha256']);source_manifest=json.loads((VACROOT/'SOURCE_MANIFEST.json').read_text());assert source_manifest['artifacts'][bindings['fmpe_estimator.pt']['path']]==a02['posterior_checkpoint_sha256']
 # Uniform spatial-core sample, declared before posterior access; all eligible rows retained.
 ix=np.sort(np.random.default_rng(20261013).choice(len(drawindex),128,replace=False));selected=[drawindex[i] for i in ix];cores=[s['core_id'] for s in selected]
 cat=fitsio.read(m['catalogue'],columns=['TARGETID','CORE_ID','Z','SUPPORTED','P_WEB']);k=np.isin(cat['CORE_ID'],cores)&cat['SUPPORTED']&(cat['Z']>=.2)&(cat['Z']<.3)
 selected_cat=cat[k];order=np.argsort(selected_cat['TARGETID']);selected_cat=selected_cat[order];ds=[];ids=[];cs=[]
 for ii,s in enumerate(selected):
  bind(s['path'],s['sha256']);mask=selected_cat['CORE_ID']==s['core_id'];ids_here=selected_cat['TARGETID'][mask]
  if not len(ids_here):continue
  with np.load(s['path']) as d:
   di=d['TARGETID'];at=np.argsort(di);loc=np.searchsorted(di[at],ids_here);assert np.array_equal(di[at[loc]],ids_here);v=d['eigenvalue_draws'][at[loc]]
  assert v.shape[1:]==(512,3)
  p=np.stack([((v>.2).sum(2)==j).mean(1) for j in range(4)],1);assert np.allclose(p,selected_cat['P_WEB'][mask],atol=1e-7,rtol=0)
  ds.append(v);ids.append(ids_here);cs.append(np.full(len(v),s['core_id']))
  if ii%32==0:print('DESI cores read',ii,flush=True)
 de=np.concatenate(ds);dids=np.concatenate(ids);dcores=np.concatenate(cs);assert len(np.unique(dids))==len(dids)==len(selected_cat)
 assert len(de)>=5000
 edges=np.linspace(-.6,1.8,121);thresholds=np.linspace(-.1,.6,71);thresholds[np.argmin(abs(thresholds-.2))]=.2
 arrays=[truth,e.reshape(-1,3),de.reshape(-1,3)];weights=[weight,np.repeat(weight/512,512),np.full(de.shape[0]*512,1/512)]
 results=[aggregate(x,w,edges,thresholds) for x,w in zip(arrays,weights)]
 hist=np.array([x[0] for x in results]);h2=np.array([x[1] for x in results]);sweep=np.array([x[2] for x in results]);at=int(np.flatnonzero(thresholds==.2)[0])
 # Exact agreement with earlier current-model mock fractions and selected DESI P_WEB.
 old=json.loads((DECK/'a07_combined/plotted_statistics.json').read_text())['class_fractions'];assert np.allclose(sweep[0,at],old['Abacus truth'],atol=1e-10);assert np.allclose(sweep[1,at],old['Abacus inferred'],atol=1e-10);assert np.allclose(sweep[2,at],selected_cat['P_WEB'].astype(float).mean(0),atol=1e-10)
 # Full support companion: same data and weights, all draw tails retained.
 limits=[float(min(x.min() for x in arrays)),float(max(x.max() for x in arrays))];full_edges=np.linspace(limits[0]-.01,limits[1]+.01,201);full=[aggregate(x,w,full_edges,thresholds) for x,w in zip(arrays,weights)]
 np.savez_compressed(out/'corner_sweep_arrays.npz',edges=edges,histograms=hist,pair_histograms=h2,thresholds=thresholds,class_fractions=sweep,full_edges=full_edges,full_histograms=np.array([x[0] for x in full]),full_pair_histograms=np.array([x[1] for x in full]),mock_evaluation_rows=rows,mock_validation_rows=vi,mock_natural_weight=weight,DESI_TARGETID=dids,DESI_CORE_ID=dcores,DESI_core_selection_index=ix)
 dump(out/'corner_sweep_provenance.json',{'model_id':a02['model_id'],'posterior_checkpoint_sha256':a02['posterior_checkpoint_sha256'],'sources':sources,'script_sha256':sha(__file__),'selection':'0.20 <= z < 0.30; supported mock evaluation rows and supported observed DESI rows. Native selections, not selection-matched cosmological populations. No CIGALE cut.','mock':{'phase':'ph006 already-exposed validation, not new blind test','rows':len(rows),'draws_per_row':512,'weights':'Saved natural_weight; each posterior draw has galaxy weight/512.'},'DESI':{'selection':'128 uniformly sampled spatial cores without replacement from all 15632 original draw-index entries; seed20261013; selection before posterior access, all eligible rows in sampled cores retained.','rows':len(de),'cores_with_eligible_rows':len(np.unique(dcores)),'draws_per_row':512,'weights':'Equal observed-galaxy weight; each posterior draw weight1/512. Same core inclusion probability for each galaxy; no inverse-selection correction.','selected_core_ids':cores},'corner':{'display_limits':[-.6,1.8],'weighted_fraction_outside_by_method_and_eigenvalue':[x[3] for x in results],'diagonal':'Weighted density histograms, 120 common bins; normalization includes out-of-window probability.','contours':'50% and 90% of total weighted pair probability, derived from common 2D histogram after Gaussian sigma=1 bin display smoothing; no fitted model.','gold_region':'Both displayed eigenvalues >0.2; not necessarily all three.','full_range_companion_limits':limits},'threshold':'Count eigenvalues strictly > threshold in every saved joint triplet; aggregate by galaxy weights. Calculated exactly via ordered-eigenvalue CDF differences, checked at0.2 against saved P_WEB.','class_fractions_at_0.2':sweep[:,at].tolist(),'full_survey_DESI_fraction_at_0.2_reference':old['DESI inferred'],'qa':{'ordering':True,'identity_join':True,'draws_reproduce_selected_DESI_P_WEB':True,'mock_fractions_reproduce_existing_figure':True,'full_range_contains_all_samples':True},'limitations':'Posterior mixtures are descriptive, not deconvolved population PDFs, volume fractions, or real-DESI calibration. DESI curve uses a spatial-core sample; differences may reflect selection, redshift mixture, sample variance, or transfer.','slurm_job':os.getenv('SLURM_JOB_ID'),'no_new_sampling_or_inference':True})
 print('Prepared corner/sweep',len(rows),len(de),'cluster fractions',sweep[:,at,3].tolist(),flush=True)

def style():
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':16,'text.color':'#F5F5F5','axes.labelcolor':'#F5F5F5','xtick.color':'#C8CDD4','ytick.color':'#C8CDD4','axes.edgecolor':'#899099','axes.facecolor':'black','figure.facecolor':'black','savefig.facecolor':'black','pdf.fonttype':42,'svg.fonttype':'none'})
def save(f,out,name,vector=True):
 for ext in (['svg','pdf','png'] if vector else ['png']):f.savefig(out/f'{name}.{ext}',dpi=240 if ext=='png' else 120)
 plt.close(f);Image.open(out/f'{name}.png').resize((1280,720),Image.Resampling.LANCZOS).save(out/f'{name}_preview.png')
def render(out):
 style();d=np.load(out/'corner_sweep_arrays.npz');r=json.loads((out/'corner_sweep_provenance.json').read_text())
 for full in [False,True]:
  edges=d['full_edges' if full else 'edges'];hist=d['full_histograms' if full else 'histograms'];h2=d['full_pair_histograms' if full else 'pair_histograms'];cent=(edges[:-1]+edges[1:])/2
  f=plt.figure(figsize=(16,9),dpi=120);gs=f.add_gridspec(3,3,left=.10,right=.66,bottom=.17,top=.83,wspace=.13,hspace=.13)
  for j in range(3):
   for i in range(j+1):
    ax=f.add_subplot(gs[j,i]);ax.set_xlim(edges[0],edges[-1]);ax.axvline(.2,color='#9CA3AD',ls=':',lw=1)
    if j==i:
     for m in range(3):ax.stairs(hist[m,i],edges,color=COLORS[m],lw=1.8,ls=STYLES[m])
     ax.set_ylim(bottom=0);ax.set_yticks([])
    else:
     pair=[(0,1),(0,2),(1,2)].index((i,j));ax.set_ylim(edges[0],edges[-1]);ax.axhline(.2,color='#9CA3AD',ls=':',lw=1);ax.fill_between([.2,edges[-1]],.2,edges[-1],color='#F2C14E',alpha=.08)
     for m in range(3):
      h=gaussian_filter(h2[m,pair],sigma=1,mode='constant');flat=np.sort(h.ravel())[::-1];cum=flat.cumsum();levels=[]
      for mass in [.9,.5]:
       if cum[-1]>=mass:levels.append(flat[min(np.searchsorted(cum,mass),len(flat)-1)])
      levels=np.unique(levels)
      if len(levels):ax.contour(cent,cent,h.T,levels=levels,colors=[COLORS[m]],linewidths=1.6,linestyles=[STYLES[m]])
    ax.spines[['top','right']].set_visible(False);ax.tick_params(labelsize=12)
    if j==2:ax.set_xlabel(['λ₁','λ₂','λ₃'][i],fontsize=21)
    else:ax.set_xticklabels([])
    if i==0 and j>0:ax.set_ylabel(['λ₁','λ₂','λ₃'][j],fontsize=21)
    elif j!=i:ax.set_yticklabels([])
  f.text(.1,.90,'Ordered tidal eigenvalues',fontsize=26);f.text(.70,.78,'Current frozen P12-A',fontsize=21)
  handles=[Line2D([],[],color=c,ls=ls,lw=2,label=n) for n,c,ls in zip(NAMES,COLORS,STYLES)];f.legend(handles=handles,loc='upper left',bbox_to_anchor=(.69,.74),frameon=False,fontsize=17,labelspacing=.7)
  f.text(.70,.50,'At λₜₕ = 0.2\nP(all three above threshold)',fontsize=17,linespacing=1.5)
  for m,c in enumerate(COLORS):f.text(.70,.40-.048*m,f"{NAMES[m]}: {100*r['class_fractions_at_0.2'][m][3]:.1f}%",fontsize=16,color=c)
  f.text(.70,.19,'Contours: 50% / 90%\nGold: both plotted λ > 0.2\nAll three methods: 0.20 ≤ z < 0.30',fontsize=14,color='#C8CDD4',linespacing=1.6)
  window='All saved-draw tails included' if full else 'Display window −0.6 to 1.8; full-range companion supplied'
  f.text(.70,.13,f"DESI: {r['DESI']['rows']:,} galaxies / {r['DESI']['cores_with_eligible_rows']} cores",fontsize=13,color='#C8CDD4');f.text(.1,.095,window,fontsize=13,color='#C8CDD4');f.text(.1,.055,'Descriptive posterior mixtures · native selections differ · real-DESI calibration unverified',fontsize=14,color='#C8CDD4')
  save(f,out,'eigenvalue_corner_full_range' if full else 'eigenvalue_corner',vector=not full)
 f,axs=plt.subplots(1,4,figsize=(16,9),dpi=120,sharey=True);f.subplots_adjust(left=.07,right=.98,bottom=.24,top=.72,wspace=.15)
 for j,ax in enumerate(axs):
  for m,ls in enumerate(STYLES):ax.plot(d['thresholds'],d['class_fractions'][m,:,j],color=WC[j],ls=ls,lw=2.2,label=NAMES[m])
  ax.axvline(.2,color='#9CA3AD',ls=':',lw=1);ax.set(xlim=(-.1,.6),ylim=(0,1),xticks=[0,.2,.4,.6]);ax.set_title(WEB[j],fontsize=24,color=WC[j],pad=18);ax.set_xlabel('Threshold λₜₕ',fontsize=18,labelpad=13);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.13)
 axs[0].set_ylabel('Expected class fraction',fontsize=20,labelpad=13)
 f.text(.5,.9,'Environment fractions depend on the threshold',ha='center',fontsize=26)
 f.text(.5,.84,'Current P12-A · 0.20 ≤ z < 0.30 · vertical line: λₜₕ = 0.2',ha='center',fontsize=17)
 handles=[Line2D([],[],color='#F5F5F5',ls=ls,lw=2,label=n) for n,ls in zip(NAMES,STYLES)];f.legend(handles=handles,loc='center',bbox_to_anchor=(.5,.79),ncol=3,frameon=False,fontsize=17)
 f.text(.5,.14,f"Abacus: {r['mock']['rows']:,} galaxies · DESI: {r['DESI']['rows']:,} galaxies in {r['DESI']['cores_with_eligible_rows']} contributing spatial cores",ha='center',fontsize=14,color='#C8CDD4')
 f.text(.5,.095,'Saved joint draws; equal weight per DESI galaxy · native sample selections differ',ha='center',fontsize=14,color='#C8CDD4');f.text(.5,.055,'Descriptive fractions, not physical volume fractions or evidence of real-DESI calibration',ha='center',fontsize=14,color='#C8CDD4')
 save(f,out,'threshold_sweep')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);p.add_argument('--stage',choices=['prepare','render','all'],default='all');a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
 if a.stage in ['prepare','all']:prepare(a.out)
 if a.stage in ['render','all']:render(a.out)
