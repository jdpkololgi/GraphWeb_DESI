#!/usr/bin/env python3
"""Matched CIGALE comparisons for measured Loa wedge environments."""
from pathlib import Path
import argparse,json,hashlib,csv
import numpy as np
from scipy.stats import rankdata,spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from desi_density_properties_20261014 import estimate,conditional_rank_corr
LABELS={'sigma5':'Projected Σ₅','n1':'Counts: 1 h⁻¹ Mpc','n2':'Counts: 2 h⁻¹ Mpc','n4':'Counts: 4 h⁻¹ Mpc','n7':'Counts: 7 h⁻¹ Mpc','trace':'Inferred trace','soft_web':'Soft web classes','concentration':'Inner / outer density','q3':'3D quadrupole','qsky':'Sky-plane quadrupole','shear':'Inferred shear proxy'}
PLABELS={'mass':'log₁₀ stellar mass [M☉]','sfr':'log₁₀ SFR [M☉ yr⁻¹]','ssfr':'log₁₀ sSFR [yr⁻¹]','quenched':'Quenched fraction','age':'log₁₀ fitted age [Gyr]','av':'Dust attenuation Aᵥ [mag]','fagn':'AGN contribution fraction'}
COLORS=['#F2C14E','#82D9CE','#A6B6D2','#C28CF5','#E85D9E'];MARKERS=['o','s','D','v','^']

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def tied_quartiles(x,strata):
 """Tie-equivalent fractional membership; every stratum has equal bin weights.
 A tied block crossing a quantile boundary has identical fractional membership
 for every tied galaxy, independent of properties or row order.
 """
 w=np.zeros((len(x),4));audit=[]
 for b in np.unique(strata):
  k=strata==b;v=x[k];lo=(rankdata(v,method='min')-1)/len(v);hi=rankdata(v,method='max')/len(v)
  a=np.maximum(0,np.minimum(hi[:,None],np.arange(1,5)[None,:]/4)-np.maximum(lo[:,None],np.arange(4)[None,:]/4))/(hi-lo)[:,None]
  assert np.allclose(a.sum(1),1) and np.allclose(a.sum(0),len(v)/4)
  w[k]=a;audit.append({'stratum':int(b),'rows':len(v),'unique_values':len(np.unique(v)),'fractional_rows':int(np.sum(np.max(a,axis=1)<1-1e-9))})
 return w,audit

def analyse(src,out):
 a=np.load(src/'raw_environment_catalogue.npz');p=json.loads((src/'measurement_provenance.json').read_text());assert sha(src/'raw_environment_catalogue.npz')==p['catalogue_sha256']
 mass=a['MASS_CG_15'];sf=a['SFR_CG_15'];age=a['AGE_CG_15'];av=a['AV_CG_15'];agn=a['FAGN_CG_15'];n=a['counts'];z=a['z'];lm=np.log10(np.maximum(mass,1e-30));ss=np.log10(np.maximum(sf,1e-30))-lm
 good=a['CIGALE_MATCHED']&(np.char.strip(a['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')&a['valid_sphere']&a['valid_sigma5']
 for x in [mass,sf,age]:good &= np.isfinite(x)&(x>0)
 good &= np.isfinite(av)&(av>=0)&np.isfinite(agn)&(agn>=0)&(agn<=1)&(lm>=8)&(lm<13)
 assert np.isfinite(a['P_WEB']).all() and np.allclose(a['P_WEB'].sum(1),1,atol=1e-6)
 eigen=a['EIGENVALUE_MEAN'].astype(float);trace=eigen.sum(1);shear=np.sqrt(.5*((eigen[:,1]-eigen[:,0])**2+(eigen[:,2]-eigen[:,1])**2+(eigen[:,2]-eigen[:,0])**2))
 vals={'sigma5':np.log10(a['sigma5']),'trace':trace,'shear':shear,'concentration':a['inner_outer_ratio'],'q3':a['Q3'],'qsky':a['Qsky']}
 # Raw counts use z-stratified quartiles; displayed x transforms preserve zeros.
 for j,rr in enumerate([1,2,4,7]):vals['n'+str(rr)]=n[:,j].astype(float)
 y={'mass':lm,'sfr':np.log10(np.maximum(sf,1e-30)),'ssfr':ss,'quenched':(ss<-11).astype(float),'age':np.log10(np.maximum(age,1e-30)/1000),'av':av,'fagn':agn}
 zb=np.floor((z-.2+1e-10)/.01).astype(int);mb=np.floor((np.nan_to_num(lm,nan=8.)-8+1e-10)/.1).astype(int);db=np.floor(np.log10(np.maximum(n[:,3],1))/.2).astype(int)
 density=['sigma5','n1','n2','n4','n7','trace','soft_web'];geometry=['concentration','q3','qsky','shear','soft_web'];shape=good&(a['Ndir3']>=20)&(a['Ndir2']>=20)&np.isfinite(a['inner_outer_ratio'])
 suites={'density':(good,density,False),'geometry':(shape,geometry,False),'geometry_at_fixed_density':(shape,geometry,True)}
 reports={};saved={'TARGETID':a['TARGETID'],'good_properties':good,'shape_selection':shape};ties={}
 for suite,(mask,methods,fixdensity) in suites.items():
  weights={};ties[suite]={};ix=np.flatnonzero(mask)
  for key in methods:
   if key=='soft_web':weights[key]=a['P_WEB'][mask].astype(float)
   else:weights[key],ties[suite][key]=tied_quartiles(vals[key][mask],zb[mask])
  stats={}
  for key,value in y.items():
   cells=zb.copy() if key=='mass' else zb+11*mb
   if fixdensity:cells=cells*100+db
   r,common=estimate(value[mask],weights,cells[mask],a['region'][mask],np.log10(np.maximum(n[mask,3],1)))
   r.pop('paired_density_minus_soft_contrast',None)
   r['control']=('redshift' if key=='mass' else 'mass + redshift')+(' + N₇ density' if fixdensity else '')
   r['label']=PLABELS[key];r['correlations']={};saved[suite+'_'+key+'_common_rows']=ix[common]
   for name in methods:
    if name=='soft_web':continue
    v=vals[name][mask];c=cells[mask];r['correlations'][name]={'spearman':float(spearmanr(v[common],value[mask][common]).statistic),'within_cell_rank':conditional_rank_corr(v[common],value[mask][common],c[common])}
    # Physical x coordinates under exactly the same common reference weights.
    xu=[]
    for g in range(4):
     by=[]
     for cell in r['common_cell_ids']:
      k=c==cell;ww=weights[name][k,g];by.append(np.sum(ww*v[k])/ww.sum())
     xu.append(float(np.dot(by,r['reference_weights'])))
    r['methods'][name]['x']=xu
   stats[key]=r
  reports[suite]={'selection_rows':int(mask.sum()),'methods':methods,'properties':stats}
  print(suite,'rows',stats['ssfr']['rows'],'cells',stats['ssfr']['cells'],'bootstrap',stats['ssfr']['bootstrap_valid'],flush=True)
  print('quenched:',{k:round(v['contrast']*100,3) for k,v in stats['quenched']['methods'].items()},flush=True)
 # Separate mass/z-only comparison on exactly the rows retained by the density control.
 # Preserve the original shape-quantile memberships rather than re-ranking this subset.
 ixshape=np.flatnonzero(shape);wshape={}
 for name in geometry:
  wshape[name]=a['P_WEB'][shape].astype(float) if name=='soft_web' else tied_quartiles(vals[name][shape],zb[shape])[0]
 same={}
 for key,value in y.items():
  ii=saved['geometry_at_fixed_density_'+key+'_common_rows'];loc=np.searchsorted(ixshape,ii)
  basecells=zb[ii] if key=='mass' else zb[ii]+11*mb[ii]
  r,common=estimate(value[ii],{n:w[loc] for n,w in wshape.items()},basecells,a['region'][ii],np.log10(n[ii,3]))
  assert common.all() and r['rows']==reports['geometry_at_fixed_density']['properties'][key]['rows']
  r.pop('paired_density_minus_soft_contrast',None);r['label']=PLABELS[key];r['control']='redshift' if key=='mass' else 'mass + redshift';same[key]=r
 reports['geometry_same_rows_without_density']={'selection_rows':'property-specific fixed-density support','methods':geometry,'properties':same}
 dump(out/'statistics.json',reports);dump(out/'tie_membership_audit.json',ties);np.savez_compressed(out/'analysis_membership.npz',**saved)
 with (out/'contrasts.csv').open('w') as f:
  w=csv.writer(f);w.writerow(['suite','property','method','contrast','interval16','interval84','rows'])
  for sn,s in reports.items():
   for pn,r in s['properties'].items():
    for mn,v in r['methods'].items():w.writerow([sn,pn,mn,v['contrast'],*v['contrast_interval'],r['rows']])
 dump(out/'analysis_provenance.json',{'measurements':p,'measurements_file':str(src/'raw_environment_catalogue.npz'),'script_sha256':sha(__file__),'estimator_dependency_sha256':sha(Path(__file__).with_name('desi_density_properties_20261014.py')),'property_rows':int(good.sum()),'shape_rows':int(shape.sum()),'controls':'Common .1dex CG15 mass/.01redshift cells and same reference mix across all estimators in each suite; mass uses z only. Geometry density control adds .2dex log10(N7) cells. Kish>=20 per group/cell. 128 paired nside8 spatial bootstrap draws, fixed strata/reference.','ties':'Four z-stratified equal-weight quantile groups. Tied counts spanning a boundary receive identical fractional membership, independently of properties; no random tie-breaking.','shear_proxy':'q=sqrt(.5*sum_pairwise(lambda_i-lambda_j)^2), evaluated at posterior-mean eigenvalues; not E[q] or a truth estimate.','Qsky_limit':'Projected directions of the same redshift-space R7 neighbours; a LOS-anisotropy sensitivity check, not a distortion-free measurement.','geometry_limit':'Q measures a quadrupole, not all angular structure. Collinear and planar configurations can both be anisotropic; clustering violates independent-neighbour assumptions in finite-N subtraction.','uncertainty':'Spatial sampling only; excludes correlated CIGALE parameter error, angular completeness, velocity distortions, and model transfer.','scope':'Same geometric wedge as preceding plots, using current uniquely matched CIGALE galaxies; not restricted to old property-cache IDs.','no_new_inference':True})

def render(out):
 R=json.loads((out/'statistics.json').read_text())
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'figure.facecolor':'#080C13','axes.facecolor':'#080C13','text.color':'#EDF1F5','axes.labelcolor':'#EDF1F5','xtick.color':'#BBC5D2','ytick.color':'#BBC5D2','axes.edgecolor':'#8995A5','pdf.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'#080C13'})
 jobs=[('local_density_comparison','density',['sigma5','n1','n2','soft_web'],'Small-scale crowding and inferred environments'),('multiscale_density_comparison','density',['n2','n4','n7','trace','soft_web'],'Density on several scales and inferred environments'),('geometry_comparison','geometry',['concentration','q3','qsky','shear','soft_web'],'Neighbour geometry and inferred environments'),('geometry_fixed_density','geometry_at_fixed_density',['concentration','q3','qsky','shear','soft_web'],'Neighbour geometry at fixed large-scale density')]
 # Also provide physical-axis versions for every measured statistic.
 for key in ['sigma5','n1','n2','n4','n7','concentration','q3','qsky']:
  suite='density' if key in ['sigma5','n1','n2','n4','n7'] else 'geometry';jobs.append((key+'_properties',suite,[key],LABELS[key]+' versus CIGALE properties'))
 for name,suite,methods,title in jobs:
  data=R[suite]['properties'];physical=len(methods)==1
  f,axs=plt.subplots(2,4,figsize=(18,10));f.subplots_adjust(left=.067,right=.978,bottom=.15,top=.78,wspace=.36,hspace=.5)
  f.suptitle(title,fontsize=23,y=.95);f.text(.5,.902,'DESI Loa sky wedge · 0.20 ≤ z < 0.30 · current 15-band CIGALE · surrounding survey neighbours',ha='center',fontsize=14)
  for ax,(key,r) in zip(axs.flat,data.items()):
   for j,m in enumerate(methods):
    v=r['methods'][m];x=np.array(v['x']) if physical else np.arange(4);mu=np.array(v['mean']);ci=np.array(v['interval']);color='#E85D9E' if m=='soft_web' else COLORS[j]
    ax.plot(x,mu,marker=MARKERS[j],lw=1.8,color=color,label=LABELS[m]);ax.vlines(x,ci[0],ci[1],color=color,lw=1.3)
   ax.set_title(PLABELS[key],fontsize=12,pad=12);ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.14)
   ax.text(.025,.965,f"{r['rows']:,} galaxies\n"+r['control'].replace(' + N₇ density','\n+ N₇ density'),transform=ax.transAxes,va='top',fontsize=9,color='#BBC5D2')
   if key=='quenched':ax.set_ylim(0,1)
   else:ax.margins(y=.35)
   if physical:
    label='Mean log₁₀ Σ₅ [h² Mpc⁻²]' if methods[0]=='sigma5' else 'Mean '+LABELS[methods[0]]
    ax.set_xlabel(label,fontsize=10)
   else:ax.set_xticks(range(4),['1','2','3','4']);ax.set_xlabel('Low → high environment group',fontsize=10)
  axs.flat[-1].axis('off');note='Identical galaxies and control mix\nfor every curve in each panel.\n\nGroups: z-stratified quantiles;\nties share fractional membership.\nSoft classes: void → wall →\nfilament → cluster.\n\nBars: 16–84% sky-block bootstrap.'
  if suite!='density':note+='\n\nShapes require ≥20 neighbours.\nFinite-N bias subtracted;\nredshift distortions remain.'
  axs.flat[-1].text(0,.96,note,va='top',fontsize=10.5,linespacing=1.5)
  if not physical:f.legend(*axs.flat[0].get_legend_handles_labels(),ncol=len(methods),loc='upper center',bbox_to_anchor=(.5,.867),frameon=False,fontsize=11)
  else:f.text(.5,.849,'Physical x values and property means use the same control weights.',ha='center',fontsize=13,color=COLORS[0])
  foot='Σ₅: fifth neighbour within ±1000 km/s · spherical radii are comoving h⁻¹ Mpc · conservative survey-edge cuts' if suite=='density' else 'Concentration: density inside 2 / density in 2–7 h⁻¹ Mpc shell · quadrupoles use neighbours inside 7 h⁻¹ Mpc'
  f.text(.5,.066,foot,ha='center',fontsize=12);f.text(.5,.027,'Same geometric wedge, newly selected current-CIGALE sample · AGE and AGN fraction are SED-fit parameters · no new ML inference',ha='center',fontsize=11,color='#BBC5D2')
  for ext in ['pdf','svg','png']:f.savefig(out/f'{name}.{ext}',dpi=180,bbox_inches='tight')
  plt.close(f);im=Image.open(out/f'{name}.png').convert('RGB');im=im.resize((1600,round(im.height*1600/im.width)),Image.Resampling.LANCZOS);im.save(out/f'{name}_preview.png')

def test():
 x=np.array([0.,0,0,0,0,0,1,2]);w,_=tied_quartiles(x,np.zeros(8));assert np.allclose(w[:6],w[0]);assert np.allclose(w.sum(0),2)
 x=np.arange(800.)%4;cells=np.repeat([0,1],400);w,_=tied_quartiles(x,cells);yy=x+3*cells;ww={'density':w,'soft_web':w}
 r,k=estimate(yy,ww,cells,np.arange(800)%13,x,minimum=10,nboot=16);assert np.isclose(r['methods']['density']['contrast'],3)
 wc,_=tied_quartiles(np.ones(800),cells);r,k=estimate(yy,{'density':wc,'soft_web':wc},cells,np.arange(800)%13,x,minimum=10,nboot=16);assert abs(r['methods']['density']['contrast'])<1e-12
 print('Tie equivalence/equal membership, preserved signal and constant-estimator null tests passed')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--source',type=Path);p.add_argument('--out',type=Path);p.add_argument('--test',action='store_true');p.add_argument('--render-only',action='store_true');a=p.parse_args()
 if a.test:test()
 else:
  a.out.mkdir(parents=True,exist_ok=True)
  if not a.render_only:analyse(a.source,a.out)
  render(a.out)
