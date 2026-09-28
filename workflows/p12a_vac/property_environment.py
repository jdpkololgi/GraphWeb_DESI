"""CIGALE wedge properties joined to unchanged current P12-A Loa posteriors."""
import argparse, hashlib, json, os
from pathlib import Path
import numpy as np
import fitsio
import healpy as hp
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm
VAC=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1/DESI_LOA_P12A_HALO48_FULL_SURVEY_VAC.fits')
CIG=Path('/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_cigale_hz/desi_wedge_env_props.parquet')
NAMES=['Void','Sheet','Filament','Knot']; COLORS=['#2585bb','#aa8310','#d25b42','#8954af']
FOOT='Provisional current-VAC environments; real-DESI coverage unverified, including this low-z wedge.\nProperty relationships are descriptive checks, not calibrated true-class population estimates or causal effects.'
def digest(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def estimate(y,p,cells,regions,nboot=96,min_eff=40):
 """Common-cell standardization and spatial cluster bootstrap at fixed reference mix."""
 if len(y)==0:return None
 nc=int(cells.max())+1;nr=int(regions.max())+1
 count=np.zeros((4,nc,nr));sums=np.zeros_like(count);sum2=np.zeros_like(count)
 flat=cells*nr+regions
 for j in range(4):
  count[j]=np.bincount(flat,weights=p[:,j],minlength=nc*nr).reshape(nc,nr)
  sums[j]=np.bincount(flat,weights=p[:,j]*y,minlength=nc*nr).reshape(nc,nr)
  sum2[j]=np.bincount(flat,weights=p[:,j]**2,minlength=nc*nr).reshape(nc,nr)
 n=count.sum(axis=2);neff=n*n/np.maximum(sum2.sum(axis=2),1e-30)
 common=np.all(neff>=min_eff,axis=0)
 if not common.any():return None
 ref=n[:,common].sum(axis=0);ref=ref/ref.sum()
 def calc(c,t):return ((t/c)*ref).sum(axis=1)
 base=calc(n[:,common],sums.sum(axis=2)[:,common]);raw=sums.sum(axis=(1,2))/count.sum(axis=(1,2))
 rng=np.random.default_rng(20260928);boots=[];rawboots=[]
 for b in range(nboot):
  mult=np.bincount(rng.integers(0,nr,nr),minlength=nr)
  c=count@mult;t=sums@mult
  if np.any(c[:,common]==0):continue
  boots.append(calc(c[:,common],t[:,common]));rawboots.append(t.sum(axis=1)/c.sum(axis=1))
 assert len(boots)>=nboot*.8
 boots=np.array(boots);rawboots=np.array(rawboots)
 keep=common[cells];w=p[keep]
 return {'mean':base.tolist(),'interval':np.quantile(boots,[.16,.84],axis=0).tolist(),'raw_mean':raw.tolist(),'raw_interval':np.quantile(rawboots,[.16,.84],axis=0).tolist(),'knot_minus_void':float(base[3]-base[0]),'contrast_interval':np.quantile(boots[:,3]-boots[:,0],[.16,.84]).tolist(),'raw_knot_minus_void':float(raw[3]-raw[0]),'rows':len(y),'common_rows':int(keep.sum()),'common_cells':int(common.sum()),'sky_blocks':nr,'effective_class_rows':((w.sum(axis=0)**2)/(w*w).sum(axis=0)).tolist(),'bootstrap_replicates':len(boots)}

def test():
 # Exactly equal within-cell properties must lose a composition-only contrast.
 cells=np.repeat([0,1],400);y=cells.astype(float)*3
 p=np.vstack([np.tile([.4,.3,.2,.1],(400,1)),np.tile([.1,.2,.3,.4],(400,1))])
 regions=np.arange(800)%13;r=estimate(y,p,cells,regions,nboot=16,min_eff=10)
 assert np.allclose(r['mean'],1.5)
 assert abs(r['raw_knot_minus_void'])>1
 r=estimate(np.ones(800)*7,p,cells,regions,nboot=16,min_eff=10)
 assert np.allclose(r['mean'],7) and np.allclose(r['raw_mean'],7)
 # With a real within-cell contrast, standardization must preserve it.
 p=np.tile(np.eye(4),(200,1));y=np.tile(np.arange(4),(200,))+cells*3
 r=estimate(y,p,cells,regions,nboot=16,min_eff=10)
 assert np.allclose(r['mean'],np.arange(4)+1.5)
 assert np.isclose(r['knot_minus_void'],3)
 print('constant-property and common-cell checks pass')

def main(out):
 out.mkdir(parents=True,exist_ok=False)
 before=digest(VAC);assert before=='cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc'
 # Explicit property whitelist: historical environment columns never enter this analysis.
 cachecols=['TARGETID','ra','dec','z','MASS_CG','SFR_CG','MASS_CG5','SFR_CG5']
 c=pd.read_parquet(CIG,columns=cachecols);assert c.TARGETID.is_unique
 # First read identities only; fetch bounded current-VAC rows after exact-ID matching.
 ids=fitsio.read(VAC,columns=['TARGETID'])['TARGETID'];order=np.argsort(ids);sortedids=ids[order]
 at=np.searchsorted(sortedids,c.TARGETID.to_numpy());matched=(at<len(ids))&(sortedids[np.minimum(at,len(ids)-1)]==c.TARGETID.to_numpy())
 cr=np.flatnonzero(matched);vr=order[at[cr]];sort=np.argsort(vr);cr=cr[sort];vr=vr[sort];c=c.iloc[cr].reset_index(drop=True)
 d=fitsio.read(VAC,rows=vr,columns=['TARGETID','RA','DEC','Z','SUPPORTED','QUALITY','P_WEB'])
 assert np.array_equal(d['TARGETID'],c.TARGETID.to_numpy())
 dz=np.abs(c.z.to_numpy()-d['Z'])/(1+d['Z']);dra=(c.ra.to_numpy()-d['RA']+180)%360-180
 sep=np.hypot(dra*np.cos(np.deg2rad(d['DEC'])),c.dec.to_numpy()-d['DEC'])*3600
 matches_before_coordinate_check=len(d)
 consistent=(dz<.001)&(sep<1)
 rejected=[{'TARGETID':int(d['TARGETID'][i]),'cached_z':float(c.z.iloc[i]),'vac_z':float(d['Z'][i]),'separation_arcsec':float(sep[i])} for i in np.flatnonzero(~consistent)]
 print('Coordinate exclusions:',rejected,flush=True)
 d=d[consistent];c=c.loc[consistent].reset_index(drop=True);dz=dz[consistent];sep=sep[consistent]
 mass=c.MASS_CG.to_numpy();sfr=c.SFR_CG.to_numpy()
 wedge=(d['RA']>=120)&(d['RA']<160)&(d['DEC']>=14.5)&(d['DEC']<30.6)&(d['Z']>=.2)&(d['Z']<.3)
 wedge_outside_ids=d['TARGETID'][~wedge].tolist()
 print('Outside current wedge:',wedge_outside_ids,flush=True)
 valid=np.isfinite(mass)&np.isfinite(sfr)&(mass>0)&(sfr>0)&d['SUPPORTED']&wedge
 lm=np.log10(mass[valid]);ss=np.log10(sfr[valid])-lm;z=d['Z'][valid];p=d['P_WEB'][valid];ra=d['RA'][valid];dec=d['DEC'][valid]
 assert np.allclose(p.sum(axis=1),1) and np.isfinite(p).all()
 assert np.all((ra>=120)&(ra<160)&(dec>=14.5)&(dec<30.6)&(z>=.2)&(z<.3))
 _,region=np.unique(hp.ang2pix(8,ra,dec,lonlat=True),return_inverse=True)
 zcell=np.floor((z-.2)/.025).astype(int)
 mb=np.floor((lm-8)/.25).astype(int);eligible=(lm>=8)&(lm<13)
 # Core standardization on the requested wedge, independently for each property.
 def run(y,weights=p,mask=eligible,controlmass=True):
  k=mask&np.isfinite(y);cs=zcell[k]+4*mb[k] if controlmass else zcell[k]
  _,cs=np.unique(cs,return_inverse=True);_,rr=np.unique(region[k],return_inverse=True)
  return estimate(y[k],weights[k],cs,rr,min_eff=20)
 results={};metrics=[('mass',lm,False,r'Mean log$_{10}$(M$_*$/M$_\odot$)'),('logssfr',ss,True,r'Mean log$_{10}$(sSFR / yr$^{-1}$)'),('low_ssfr',(ss<-11).astype(float),True,'Low-sSFR fraction (< −11 dex)')]
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.15,'figure.facecolor':'white','savefig.facecolor':'white'})
 pdf=PdfPages(out/'LOA_CIGALE_environment.pdf');files=[]
 def save(fig,name):
  fig.text(.5,.012,FOOT,ha='center',va='bottom',fontsize=9,color='#555555');fig.savefig(out/(name+'.png'),dpi=180);pdf.savefig(fig);plt.close(fig);files.append(name+'.png')
 def draw(a,r):
  a.plot(range(4),r['raw_mean'],'--',color='#999999',label='Raw mean')
  a.plot(range(4),r['mean'],color='#333333',lw=1,label='Controlled mean')
  low,high=np.array(r['interval'])
  for j in range(4):
   a.vlines(j,low[j],high[j],color=COLORS[j],lw=2);a.scatter(j,r['mean'][j],color=COLORS[j],s=45,zorder=4)
  a.set_xticks(range(4),NAMES);a.set_xlim(-.3,3.3)
 fig,axs=plt.subplots(1,3,figsize=(15,6.3));fig.subplots_adjust(left=.07,right=.98,top=.74,bottom=.23,wspace=.3)
 fig.suptitle('Galaxy properties versus current Loa environment',fontsize=21,y=.97)
 fig.text(.5,.85,'Saved CIGALE CG15 properties · RA 120–160°, Dec 14.5–30.6°, 0.20 ≤ z < 0.30\nGrey dashed: raw probability-weighted mean. Coloured: common z mix (mass); common z + CIGALE mass mix (sSFR / fraction).\nIntervals: 16–84% sky-block bootstrap; property errors and model-transfer systematics are not propagated.',ha='center',fontsize=10)
 for a,(name,y,mc,label) in zip(axs,metrics):
  r=run(y,controlmass=mc);assert r;results[name]=r;draw(a,r);a.set_ylabel(label);a.set_title(f"N={r['common_rows']:,} retained; {r['common_cells']} strata",fontsize=11)
 axs[0].legend(loc='best',fontsize=9)
 save(fig,'01_cigale_environment_trends')
 # Fixed mass-bin curves with redshift balancing within each bin.
 fig,axs=plt.subplots(1,2,figsize=(13,6.4));fig.subplots_adjust(left=.08,right=.98,top=.79,bottom=.24,wspace=.3)
 fig.suptitle('Does the relationship survive at fixed stellar mass?',fontsize=20,y=.97)
 fig.text(.5,.865,'Current environment probabilities · each mass bin shares a common redshift distribution\nShading: 16–84% sky-block bootstrap; bins without common support are omitted.',ha='center',fontsize=11)
 edges=np.arange(9,12.51,.5);centres=(edges[:-1]+edges[1:])/2
 for a,(name,y,label) in zip(axs,[('logssfr',ss,r'Mean log$_{10}$(sSFR / yr$^{-1}$)'),('low_ssfr',(ss<-11).astype(float),'Low-sSFR fraction (< −11 dex)')]):
  rr=[run(y,mask=eligible&(lm>=lo)&(lm<hi),controlmass=False) for lo,hi in zip(edges[:-1],edges[1:])];results[name+'_by_mass']=rr
  for j in range(4):
   mean=np.array([r['mean'][j] if r else np.nan for r in rr]);low=np.array([r['interval'][0][j] if r else np.nan for r in rr]);high=np.array([r['interval'][1][j] if r else np.nan for r in rr])
   a.plot(centres,mean,color=COLORS[j],marker=['o','s','^','D'][j],label=NAMES[j]);a.fill_between(centres,low,high,color=COLORS[j],alpha=.13)
  a.set(xlabel=r'CIGALE log$_{10}$(M$_*$/M$_\odot$)',ylabel=label);a.legend(fontsize=10)
 save(fig,'02_fixed_mass_relations')
 # Weighted bivariate distributions; identical histogram ranges and normalization across environments.
 fig,axs=plt.subplots(2,2,figsize=(12,10));fig.subplots_adjust(left=.09,right=.86,top=.85,bottom=.17,hspace=.25,wspace=.2)
 fig.suptitle('Stellar mass and specific star formation',fontsize=21,y=.97)
 fig.text(.5,.905,'Same galaxies in each panel, weighted by their current environment probabilities\nDensity normalized by total class probability weight; identical colour scales. These are observed selected-sample distributions.',ha='center',fontsize=10)
 xe=np.linspace(8.5,12.5,65);ye=np.linspace(-14,-8,81);histstats=[]
 for j,a in enumerate(axs.flat):
  hist,_,_=np.histogram2d(lm,ss,bins=[xe,ye],weights=p[:,j]);den=p[:,j].sum();density=hist/(den*np.diff(xe)[:,None]*np.diff(ye)[None,:]);im=a.pcolormesh(xe,ye,density.T,cmap='magma',norm=LogNorm(vmin=.001,vmax=2),rasterized=True)
  a.axhline(-11,color='#43cbd5',ls='--',lw=1);a.set(title=NAMES[j],xlabel=r'log$_{10}$(M$_*$/M$_\odot$)',ylabel=r'log$_{10}$(sSFR / yr$^{-1}$)');histstats.append({'class':NAMES[j],'class_weight':float(den),'in_plot_fraction':float(hist.sum()/den)})
 cax=fig.add_axes([.89,.25,.018,.51]);fig.colorbar(im,cax=cax,label='Probability-weighted density [dex⁻²]')
 save(fig,'03_mass_ssfr_distributions')
 # Distinguish uncertain-membership dilution and CIGALE configuration sensitivity.
 hard=np.eye(4)[p.argmax(axis=1)];high=p.max(axis=1)>=.8
 c5m=c.MASS_CG5.to_numpy()[valid];c5s=c.SFR_CG5.to_numpy()[valid];good5=(c5m>0)&(c5s>0)&np.isfinite(c5m)&np.isfinite(c5s)
 ss5=np.full(len(lm),np.nan);ss5[good5]=np.log10(c5s[good5]/c5m[good5])
 # CG5 control uses its own masses in matching cells.
 def config5(y):
  k=good5;mass5=np.log10(c5m[k]);kix=np.flatnonzero(k)[(mass5>=8)&(mass5<13)];mass5=np.log10(c5m[kix]);cs=zcell[kix]+4*np.floor((mass5-8)/.25).astype(int)
  _,cs=np.unique(cs,return_inverse=True);_,rr=np.unique(region[kix],return_inverse=True)
  return estimate(y[kix],p[kix],cs,rr,min_eff=20)
 fig,axs=plt.subplots(1,2,figsize=(12,6));fig.subplots_adjust(left=.24,right=.98,top=.77,bottom=.24,wspace=.4)
 fig.suptitle('Knot − void contrast: membership and property checks',fontsize=18,y=.97)
 fig.text(.5,.865,'All comparisons control redshift and stellar mass; retained supports differ by estimator\nIntervals: 16–84% spatial bootstrap. High-confidence selection is a sensitivity check, not the primary sample.',ha='center',fontsize=10)
 for a,(name,y,_,label) in zip(axs,metrics[1:]):
  rrs=[results[name],run(y,weights=hard),run(y,weights=hard,mask=eligible&high),config5(ss5 if name=='logssfr' else (ss5<-11).astype(float))]
  results[name+'_sensitivity']=rrs
  for i,r in enumerate(rrs):
   if r:a.hlines(i,*r['contrast_interval'],color='#2585bb',lw=2);a.plot(r['knot_minus_void'],i,'o',color='#2585bb')
  a.set_yticks(range(4),['CG15 probability weights','CG15 argmax','CG15 argmax; max P ≥ 0.8','CG5 probability weights']);a.invert_yaxis();a.axvline(0,color='grey',ls=':');a.set_xlabel('Knot − void: '+('mean log sSFR [dex]' if name=='logssfr' else 'low-sSFR fraction'))
 axs[1].set_yticklabels([])
 save(fig,'04_membership_sensitivity');pdf.close()
 # Audit availability by inferred class on every matched supported wedge object.
 allp=d['P_WEB'][d['SUPPORTED']&wedge];validp=d['P_WEB'][valid]
 coverage=(validp.sum(axis=0,dtype='f8')/allp.sum(axis=0,dtype='f8')).tolist()
 report={'vac':str(VAC),'vac_sha256':before,'vac_unchanged':digest(VAC)==before,'cache':str(CIG),'cache_sha256':digest(CIG),'cache_columns_read':cachecols,'cache_original_rows':int(matched.size),'current_vac_matches_before_coordinate_check':matches_before_coordinate_check,'coordinate_exclusions':rejected,'current_vac_matches':len(d),'supported_matches':int(d['SUPPORTED'].sum()),'outside_current_wedge_targetids':wedge_outside_ids,'supported_current_wedge_matches':int((d['SUPPORTED']&wedge).sum()),'positive_finite_mass_sfr_rows':int(valid.sum()),'mass_range_8_13_rows':int(eligible.sum()),'max_normalized_redshift_difference':float(dz.max()),'max_position_difference_arcsec':float(sep.max()),'class_property_availability':coverage,'class_effective_rows':((p.sum(axis=0)**2)/(p*p).sum(axis=0)).tolist(),'sky_blocks_nside8':int(region.max()+1),'results':results,'histogram_coverage':histstats,'figures':files,'method':'Descriptive probability-weighted means, not true-class population estimators. z bins .025; CIGALE mass bins .25 dex over 8<=logM<13; common cells require >=20 effective rows in every class; reference is pooled retained sample. Raw uses full eligible sample. 96 nside8 block resamples, common reference mix held fixed; percentile16-84. No calibrated spatial posterior, property-error or transfer-systematic propagation.','selection':'Original VAC unchanged. Analysis requires cached CIGALE goodPhoto match, finite positive SFR and mass; nonpositive SFR not assigned quenched by default. Low-sSFR threshold -11 is a descriptive proxy.','script_sha256':digest(__file__),'slurm_job_id':os.environ.get('SLURM_JOB_ID'),'science_release_ready':False}
 assert report['vac_unchanged']
 (out/'RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'ARTIFACT_SHA256.json').write_text(json.dumps({f.name:digest(f) for f in out.iterdir() if f.is_file()},indent=2)+'\n')
 print(json.dumps({k:report[k] for k in ['current_vac_matches','positive_finite_mass_sfr_rows','class_property_availability','sky_blocks_nside8']},indent=2));print(json.dumps({n:results[n] for n in ['mass','logssfr','low_ssfr']},indent=2),flush=True)

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--test',action='store_true');ap.add_argument('--out',type=Path);a=ap.parse_args()
 if a.test:test()
 else:main(a.out)
