"""Bounded-column Loa census; empirical BRIGHT completeness, no model inference."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
import numpy as np, fitsio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from legacy_selection_tests_20261014 import completeness
H=Path('/global/u2/d/dkololgi');R=H/'fmsu_tokyo_20261010';O=R/'loa_mass_redshift_completeness_20261015';O.mkdir(exist_ok=True)
VAC=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')
PHOT=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
COLS=['TARGETID','RA','DEC','Z','SUPPORTED','CIGALE_MATCHED','CIGALE_AMBIGUOUS','CIGALE_SPECTYPE','MASS_CG_15']
def sky(a):return (a['RA']>=120)&(a['RA']<160)&(a['DEC']>=14.5)&(a['DEC']<30.6)
def save(fig,name):
 fig.savefig(O/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(O/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
parts=[];census={}
with fitsio.FITS(VAC) as f:
 census['vac_rows']=f[1].get_nrows()
 for start in range(0,census['vac_rows'],200000):
  d=f[1].read(rows=np.arange(start,min(start+200000,census['vac_rows'])),columns=COLS)
  good=d['SUPPORTED'].astype(bool)&d['CIGALE_MATCHED'].astype(bool)&~d['CIGALE_AMBIGUOUS'].astype(bool)&(np.char.strip(d['CIGALE_SPECTYPE'].astype(str))=='GALAXY')&np.isfinite(d['MASS_CG_15'])&(d['MASS_CG_15']>0)&(d['Z']>=.15)&(d['Z']<.55)
  parts.append(d[good])
a=np.concatenate(parts);del parts;print('Selected VAC rows',len(a),flush=True)
order=np.argsort(a['TARGETID']);sid=a['TARGETID'][order];assert np.all(np.diff(sid)>0)
r=np.full(len(a),np.nan);lim=r.copy();found=np.zeros(len(a),np.uint8);sys=np.full(len(a),'',dtype='U1');checks=0
with fitsio.FITS(PHOT) as f:
 census['photometry_rows']=f[1].get_nrows()
 for start in range(0,census['photometry_rows'],200000):
  d=f[1].read(rows=np.arange(start,min(start+200000,census['photometry_rows'])),columns=['TARGETID','RA','DEC','Z_not4clus','FLUX_R','MW_TRANSMISSION_R','PHOTSYS'])
  at=np.searchsorted(sid,d['TARGETID']);ok=(at<len(sid))&(sid[np.minimum(at,len(sid)-1)]==d['TARGETID']);ii=order[at[ok]];dd=d[ok]
  for name,other in [('RA','RA'),('DEC','DEC'),('Z_not4clus','Z')]:assert np.array_equal(dd[name],a[other][ii]);checks+=1
  flux=dd['FLUX_R']/dd['MW_TRANSMISSION_R'];ps=np.char.strip(dd['PHOTSYS'].astype(str));assert np.isin(ps,['N','S']).all();good=np.isfinite(flux)&(flux>0)
  r[ii[good]]=22.5-2.5*np.log10(flux[good]);lim[ii]=np.where(ps=='N',19.54,19.5);sys[ii]=ps;found[ii]+=1
assert np.all(found==1);print('Photometry join complete',flush=True)
mass=np.log10(a['MASS_CG_15']);bright=np.isfinite(r)&(r<lim);census['supported_valid_mass_rows']=len(a);census['explicit_bright_rows']=int(bright.sum());census['invalid_or_outside_bright']=int((~bright).sum())
old=np.load(R/'legacy_gat_comparison_20261014/legacy_cigale.npz');om=np.log10(np.maximum(old['MASS_CG_15'],1e-30));og=np.isfinite(om)&(old['MASS_CG_15']>0)&(np.char.strip(old['SPECTYPE'].astype(str))=='GALAXY')
# Original low-z parent is needed BEFORE the old fixed mass cut.
OP=Path('/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/loa-combined-lowz.fits')
parent=fitsio.read(OP,columns=['TARGETID','RA','DEC','Z','ZWARN','DELTACHI2','BGS_TARGET','SURVEY','FLUX_R','LOGMSTAR'])
_,inv,n=np.unique(parent['TARGETID'],return_inverse=True,return_counts=True)
op=(n[inv]==1)&(parent['ZWARN']==0)&(parent['DELTACHI2']>=25)&((parent['BGS_TARGET']&2)!=0)&(np.char.strip(parent['SURVEY'].astype(str))=='main')&(parent['FLUX_R']>0)&np.isfinite(parent['LOGMSTAR'])
oldmass=parent['LOGMSTAR']-2*np.log10(.704);oldr=22.5-2.5*np.log10(np.maximum(parent['FLUX_R'],1e-30));op&=oldr<19.5
results={'census':census,'scopes':{},'method':'Faintest20% in r-r_limit in each dz=.01 bin;90th/95th percentile logM+.4(r-r_limit). No mass cutoff before estimation; no extrapolation;minimum100 objects/bin. Old FSF parent dz=.005 and historical flux convention. Current CG15 masses;BRIGHT limitsN19.54/S19.5 with extinction correction.','source_paths':[str(VAC),str(PHOT),str(OP)],'source_stats':{str(p):{'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in [VAC,PHOT,OP]},'checks':{'unique_current_ids':True,'all_photometry_matches_once':True,'exact_id_ra_dec_z_parity':True,'chunk_coordinate_checks':checks},'limits':'Conditional on valid CIGALE matches and SUPPORTED predictions, not total photometric/spectroscopic completeness. Old saved CG15 rows already mass-selected in FSF. No current values below.15 and no interpolation across gap. Curves estimated, not calibrated guarantees.'}
for scope in ['full_footprint','property_wedge']:
 use=bright if scope=='full_footprint' else bright&sky(a);ou=og if scope=='full_footprint' else og&sky(old);oldparent=op if scope=='full_footprint' else op&sky(parent)
 curve=completeness(a['Z'][use],mass[use],r[use],lim[use],np.round(np.arange(.15,.55001,.01),8));oc=completeness(parent['Z'][oldparent],oldmass[oldparent],oldr[oldparent],np.full(oldparent.sum(),19.5),np.linspace(.01,.06000001,11))
 h,xe,ye=np.histogram2d(a['Z'][use],mass[use],bins=[np.linspace(0,.56,281),np.linspace(7,13,241)]);ho,_,_=np.histogram2d(old['z'][ou],om[ou],bins=[xe,ye]);hist=h+ho
 np.savez_compressed(O/(scope+'_histogram.npz'),current=h,old_saved=ho,z_edges=xe,mass_edges=ye)
 stats=[]
 for c in curve:
  sel=use&(a['Z']>=c['zlo'])&(a['Z']<c['zhi']);c['fraction_above_1e9']=float(np.mean(mass[sel]>=9)) if sel.any() else None;c['mass_median']=float(np.median(mass[sel])) if sel.any() else None
  if c['usable']:c['fraction_above_limit']=float(np.mean(mass[sel]>=c['logmass90']))
 results['scopes'][scope]={'current_rows':int(use.sum()),'old_saved_cg15_rows':int(ou.sum()),'current_curve':curve,'old_fsf_curve':oc,'histogram_in_range_rows':int(hist.sum())}
 plt.rcParams.update({'font.size':11,'axes.spines.top':True,'axes.spines.right':True})
 fig,ax=plt.subplots(figsize=(11.5,7));mesh=ax.pcolormesh(xe,ye,np.ma.masked_equal(hist.T,0),norm=LogNorm(vmin=1,vmax=max(2,hist.max())),cmap='magma_r',rasterized=True)
 ax.axvspan(.06,.15,color='#d2d7dc',alpha=.65,zorder=2);ax.axvspan(.01,.06,color='#24a0ad',alpha=.08)
 ax.axvline(.06,color='#188e9b',lw=2);ax.axvline(.15,color='#656565',ls='--',lw=1.8);ax.axvline(.55,color='#656565',ls=':',lw=1.5)
 good=[c for c in curve if c['usable']];zc=[(c['zlo']+c['zhi'])/2 for c in good]
 ax.plot(zc,[c['logmass90'] for c in good],color='#e600ad',lw=2.5,label='BGS Bright: estimated 90% mass limit')
 ax.plot(zc,[c['logmass95'] for c in good],color='#9b55ad',lw=1.8,ls='--',label='95th-percentile sensitivity curve')
 ax.axhline(9,color='#138b95',ls='-.',lw=1.8,label='10⁹ M☉ reference — not completeness at high z')
 ax.text(.032,12.86,'Old GAT\n0.01–0.06',ha='center',va='top',fontsize=10,color='#086770',bbox=dict(fc='white',alpha=.9,ec='none'))
 ax.text(.105,10.7,'No overlap\nbetween\nsaved VACs',ha='center',va='center',fontsize=9,color='#444444')
 ax.annotate('P12-A begins at z = 0.15\n(low-z mock truth-mapping restriction)',xy=(.15,8.15),xytext=(.23,7.4),arrowprops=dict(arrowstyle='->',color='#555555'),fontsize=10)
 ax.axvspan(.2,.3,ymin=0,ymax=.025,color='#ef9927');ax.text(.25,7.22,'Earlier property tests',ha='center',fontsize=9)
 ax.set(xlim=(0,.56),ylim=(7,13),xlabel='Redshift z',ylabel='log₁₀(stellar mass / M☉), CIGALE CG15',title='Loa mass–redshift distribution: old GAT versus current P12-A\n'+('Full available footprint' if scope=='full_footprint' else 'Property-comparison wedge: RA 120–160°, Dec 14.5–30.6°'))
 ax.legend(loc='upper left',bbox_to_anchor=(.32,.99),fontsize=9,framealpha=.95);fig.colorbar(mesh,ax=ax,pad=.02,label='Galaxies per mass–redshift bin (log colour scale)')
 ax.tick_params(which='both',direction='in',top=True,right=True);ax.minorticks_on()
 fig.text(.5,-.085,f'Current: {use.sum():,} SUPPORTED, uniquely CIGALE-matched Bright galaxies; old: {ou.sum():,} matched saved-GAT galaxies.\nOld masses shown on the same CG15 scale, but its historical selection used FastSpecFit: no exact CG15 mass box is implied.\nCurves use the faintest 20% before any stellar-mass cut. No extrapolation below z=0.15 or beyond z=0.55.',ha='center',fontsize=9)
 fig.tight_layout();save(fig,'loa_mass_z_'+scope)
 # Separate original-convention panel to avoid drawing an FSF threshold on CG15 masses.
 fig,axs=plt.subplots(1,2,figsize=(12,5.3))
 hh,xx,yy=np.histogram2d(parent['Z'][oldparent],oldmass[oldparent],bins=[np.linspace(.01,.061,70),np.linspace(7,12.5,150)])
 axs[0].pcolormesh(xx,yy,np.ma.masked_equal(hh.T,0),norm=LogNorm(vmin=1,vmax=max(2,hh.max())),cmap='magma_r',rasterized=True)
 cc=[c for c in oc if c['usable']];axs[0].plot([(c['zlo']+c['zhi'])/2 for c in cc],[c['logmass90'] for c in cc],color='#e600ad',lw=2,label='Reconstructed 90% limit')
 threshold=9-2*np.log10(.704);axs[0].axhline(threshold,color='#1597a4',lw=2,label=f'Old raw FSF ≥9 → physical ≥{threshold:.3f}');axs[0].axhline(9,c='gray',ls=':',label='Intended physical 10⁹ reference');axs[0].axvline(.06,color='#1597a4',lw=2)
 axs[0].set(xlim=(.009,.063),ylim=(7.5,12.5),xlabel='Redshift z',ylabel='log₁₀ stellar mass: FastSpecFit, h=0.704',title='Old selection: fixed mass cut + nearby ceiling');axs[0].legend(fontsize=8,loc='upper left')
 axs[1].pcolormesh(xe,ye,np.ma.masked_equal(h.T,0),norm=LogNorm(vmin=1,vmax=max(2,h.max())),cmap='magma_r',rasterized=True)
 axs[1].plot(zc,[c['logmass90'] for c in good],c='#e600ad',lw=2,label='Current estimated 90% limit');axs[1].axhline(9,c='#1597a4',ls='-.',label='10⁹ reference (not complete)');axs[1].set(xlim=(.15,.55),ylim=(7.5,12.5),xlabel='Redshift z',ylabel='log₁₀ stellar mass: CIGALE CG15',title='Current selection: brightness-limited, z=0.15–0.55');axs[1].legend(fontsize=8,loc='lower right')
 fig.suptitle('Why the two Loa selections describe different populations — '+scope.replace('_',' '));fig.tight_layout(rect=[0,.05,1,.94]);fig.text(.5,.005,'Different mass estimators are explicitly separated. Old completeness uses the pre-mass-cut main-BRIGHT parent.\nThe z=0.15 boundary reflects mock-label validity, not a stellar-mass completeness threshold.',ha='center',fontsize=9);save(fig,'selection_conventions_'+scope)
 print(scope,'done',use.sum(),flush=True)
(O/'RESULTS.json').write_text(json.dumps(results,indent=2,allow_nan=False));(O/'MANIFEST.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file() and p.name!='MANIFEST.json'},indent=2))
print('COMPLETE',flush=True)
