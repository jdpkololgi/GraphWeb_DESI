"""Filter-integrated smooth continuum and profile-proxy BGS selection tests."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import sys,json,hashlib
from pathlib import Path
import numpy as np
import h5py
from scipy.interpolate import interp1d
from scipy.special import gammaincinv
from scipy.stats import ncx2
from astropy.cosmology import Planck18
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'/global/common/software/desi/perlmutter/desiconda/current/code/speclite/v1.0.0/lib/python3.13/site-packages')
from speclite import filters
from desitarget.cuts import isBGS_colors
from tng_bright_selection_20261015 import quarters
from loa_tng_abacus_overlay_20261015 import property_data
R=Path('/global/u2/d/dkololgi/fmsu_tokyo_20261010');T=R/'tng_density_control_20261015';O=R/'tng_bright_refined_20261015';O.mkdir(exist_ok=True)
S=Path('/pscratch/sd/d/dkololgi/TNG300-1/groups_099')

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def photometry():
 old=np.load(R/'tng_bright_selection_20261015/photometry.npz');sid=old['subhalo_id'];cache=O/'inputs.npz'
 if not cache.exists():
  mags=np.full((len(sid),4),np.nan);size=np.full(len(sid),np.nan);offset=0
  for j in range(600):
   with h5py.File(S/f'fof_subhalo_tab_099.{j}.hdf5') as f:
    n=int(f['Header'].attrs['Nsubgroups_ThisFile']);ii=np.flatnonzero((sid>=offset)&(sid<offset+n))
    if len(ii):
     local=sid[ii]-offset;mags[ii]=f['Subhalo/SubhaloStellarPhotometrics'][:,4:8][local]
     size[ii]=f['Subhalo/SubhaloHalfmassRadType'][:,4][local]/float(f['Header'].attrs['HubbleParam'])
    offset+=n
   if j%100==0:print('input chunk',j,flush=True)
  assert np.isfinite(mags).all() and np.all(size>0)
  assert np.allclose(mags[:,:2],np.array([old['g'],old['r']]).T)
  np.savez_compressed(cache,ids=old['ids'],subhalo_id=sid,mags=mags,rhalf_pkpc=size)
 return np.load(cache)

def filter_weights(names,wave):
 fs=filters.load_filters(*names);out=[]
 for f in fs:
  w=np.interp(wave,f.wavelength,f.response,left=0,right=0)*np.gradient(wave)/wave
  out.append(w/w.sum())
 return np.array(out),np.array([f.effective_wavelength.value for f in fs])

def build_observed(p):
 cache=O/'observed_photometry.npz'
 if cache.exists():return np.load(cache),json.loads((O/'photometry_validation.json').read_text())
 wave=np.arange(2900.,11501.,10.)
 sw,knots=filter_weights(['sdss2010-'+x for x in 'griz'],wave)
 dw,_=filter_weights(['decam2014-'+x for x in 'grz'],wave)
 def basis(w):return interp1d(np.log(knots),np.eye(4),axis=0,fill_value='extrapolate')(np.log(w)).T
 b0=basis(wave);n=len(p['ids']);out=np.empty((10,n,3));err=[];nodelogs=np.empty((n,4))
 # Fit a positive piecewise power-law f_nu continuum. The nodal fluxes are
 # iteratively adjusted so integrated SDSS photometry matches all four inputs.
 for start in range(0,n,1024):
  end=min(start+1024,n);target=-.4*p['mags'][start:end];nodes=target.copy()
  for iteration in range(25):
   pred=np.log10(10**(nodes@b0)@sw.T);nodes+=.8*(target-pred)
  err.extend(np.abs(2.5*(pred-target)).ravel());nodelogs[start:end]=nodes
  for j in range(10):
   z=.205+.01*j;flux=10**(nodes@basis(wave/(1+z)))@dw.T
   out[j,start:end]=-2.5*np.log10(flux*(1+z))+Planck18.distmod(z).value
 print('continuum reconstructed',np.quantile(err,[.5,.95,.99,1]),flush=True)
 assert max(err)<.025,'Continuum failed input photometry reconstruction'
 # A flat f_nu spectrum must give K=-2.5log10(1+z) for any filter.
 assert np.allclose((np.ones(len(wave))@dw.T),1)
 # Fiber aperture: circular Sersic with Re=.75*3D stellar half-mass radius,
 # 1-arcsec FWHM Gaussian seeing, 1.5-arcsec diameter fiber. Neither the
 # projected half-light radius nor Sersic index is measured by this catalogue.
 reg=np.geomspace(.005,100,350);sigma=1/2.354820045;frac={}
 prob=(np.arange(600)+.5)/600
 for sn in [1,4]:
  radial=(gammaincinv(2*sn,prob)/gammaincinv(2*sn,.5))**sn
  fr=ncx2.cdf((.75/sigma)**2,2,(reg[:,None]*radial[None,:]/sigma)**2).mean(1)
  assert np.all((fr>0)&(fr<1)) and np.all(np.diff(fr)<=1e-10)
  frac[sn]=np.array([np.interp(.75*p['rhalf_pkpc']/Planck18.kpc_proper_per_arcmin(.205+.01*j).value*60,reg,fr) for j in range(10)])
 np.savez_compressed(cache,mag_grz=out,fiber_n1=frac[1],fiber_n4=frac[4],logfnu_nodes=nodelogs)
 val={'input_reconstruction_abs_mag_error_quantiles':np.quantile(err,[.5,.95,.99,1]).tolist(),'flat_fnu_normalization':True,'fiber_fraction_monotonic':True,'method':'Four-band constrained positive smooth f_nu continuum; integrated SDSS2010 input and DECam2014 output response curves, 10A quadrature. Not a stellar-population template or dust model. Sersic n1/n4 profile proxy with Re=.75 Rhalf3D; 1arcsec Gaussian seeing and .75arcsec aperture radius.'}
 (O/'photometry_validation.json').write_text(json.dumps(val,indent=2));return np.load(cache),val

def main():
 p=photometry();obs,pval=build_observed(p);ids=p['ids'];d=np.load(T/'catalogue.npz');m=np.log10(d['mass_msun'][ids]);q=(d['sfr'][ids]/d['mass_msun'][ids]<1e-11).astype(float)
 tr=np.load(T/'truth.npz')['eig_R7_n512'][ids].sum(1);block=np.ravel_multi_index(np.minimum((d['pos_hmpc'][ids]/205*5).astype(int),4).T,(5,5,5))
 a=np.load(R/'loa_tng_abacus_overlay_20261015/loa_predictions.npz');ph=np.load(R/'legacy_selection_tests_20261014/p12a_targeting.npz');assert np.array_equal(a['TARGETID'],ph['TARGETID'])
 good,y=property_data(a);good&=(ph['r']<ph['rlimit'])&(y['mass']<11.5)&(a['z']>=.2)&(a['z']<.3)
 lm=y['mass'][good];lz=np.floor((a['z'][good]-.2)/.01+1e-8).astype(int);lq=y['quenched'][good];lt=a['abacus'][good].sum(1);_,lb=np.unique(a['region'][good],return_inverse=True)
 td=np.random.default_rng(731).multinomial(125,np.ones(125)/125,200);nb=lb.max()+1;ld=np.random.default_rng(732).multinomial(nb,np.ones(nb)/nb,200)
 arms=['uncut','old_powerlaw','filter_r','bgs_n1','bgs_n4'];masks={};counts=[]
 for j in range(10):
  z=.205+.01*j;oldr=p['mags'][:,1]+Planck18.distmod(z).value+(p['mags'][:,0]-p['mags'][:,1])*np.log(1+z)/np.log(6231/4770)-2.5*np.log10(1+z)
  mags=obs['mag_grz'][j];flux=10**((22.5-mags)/2.5)
  masks[j]={'uncut':np.ones(len(ids),bool),'old_powerlaw':oldr<19.5,'filter_r':mags[:,1]<19.5}
  for n in [1,4]:
   ff=flux[:,1]*obs[f'fiber_n{n}'][j]
   mask=isBGS_colors(gflux=flux[:,0],rflux=flux[:,1],zflux=flux[:,2],w1flux=np.ones(len(ids)),rfiberflux=ff,rfibertotflux=ff,south=True,targtype='bright')
   # Independent expression of applicable southern bright colour/fiber cuts.
   rf=mags[:,1]-2.5*np.log10(obs[f'fiber_n{n}'][j]);g,r,zz=mags.T
   direct=(g-r>-1)&(g-r<4)&(r-zz>-1)&(r-zz<4)&(r<19.5)&(r>=12)&(rf>=15)&(((r<17.8)&(rf<r+5.1))|((r>17.8)&(r<20)&(rf<22.9)))
   assert np.array_equal(mask,direct)
   assert np.all(~mask|masks[j]['filter_r']);masks[j][f'bgs_n{n}']=mask
  counts.append({'z':z,**{k:int(v.sum()) for k,v in masks[j].items()},'r_refined_minus_old_quantiles':np.quantile(mags[:,1]-oldr,[.1,.5,.9]).tolist()})
 checks=0
 def stat(yy,xx,bb,draw):
  nonlocal checks
  w=quarters(xx);nb=draw.shape[1];den=np.array([np.bincount(bb,weights=w[:,k],minlength=nb) for k in range(4)]).T;num=np.array([np.bincount(bb,weights=w[:,k]*yy,minlength=nb) for k in range(4)]).T
  v=num.sum(0)/den.sum(0);assert np.allclose(v.mean(),yy.mean());checks+=1
  db=draw@den; boot=np.divide(draw@num,db,out=np.full_like(db,np.nan),where=db>0)
  return v,boot
 results={};rows=[]
 fig,axs=plt.subplots(2,2,figsize=(12,9),sharex=True,sharey=True)
 for row,width in enumerate([.25,.1]):
  mb=np.floor((m-9)/width+1e-9).astype(int);lmb=np.floor((lm-9)/width+1e-9).astype(int);stats={k:{} for k in arms+['loa']};cells=[]
  for j in range(10):
   for b in range(round(2.5/width)):
    li=np.flatnonzero((lmb==b)&(lz==j));ti={k:np.flatnonzero((mb==b)&mask) for k,mask in masks[j].items()}
    if len(li)<80 or min(map(len,ti.values()))<80:continue
    key=(b,j);cells.append({'key':key,'loa':len(li),'tng_counts':{k:len(ii) for k,ii in ti.items()}})
    stats['loa'][key]=stat(lq[li],lt[li],lb[li],ld)
    for k,ii in ti.items():stats[k][key]=stat(q[ii],tr[ii],block[ii],td)
  for col,high in enumerate([False,True]):
   selected=[c for c in cells if not high or 9+c['key'][0]*width>=10.5-1e-8];weights=np.array([c['loa'] for c in selected],float);weights/=weights.sum();curves={};boots={}
   for arm in arms+['loa']:
    v=sum(w*stats[arm][tuple(c['key'])][0] for c,w in zip(selected,weights));bt=sum(w*stats[arm][tuple(c['key'])][1] for c,w in zip(selected,weights));boots[arm]=bt
    assert np.all(np.isfinite(bt).sum(0)>=160), 'Too few valid spatial replicates'
    curves[arm]={'valid_bootstrap_min':int(np.isfinite(bt).sum(0).min()),'mean':v.tolist(),'ci95':np.nanquantile(bt,[.025,.975],axis=0).tolist(),'contrast_pp':float(100*(v[-1]-v[0]))}
   for arm in arms:
    gap=100*(np.array(curves[arm]['mean'])-curves['loa']['mean']);curves[arm]['mean_gap_pp']=float(gap.mean());curves[arm]['mean_gap_ci95_pp']=(100*np.nanquantile((boots[arm]-boots['loa']).mean(1),[.025,.975])).tolist()
   name=f'{width}_{"high" if high else "full"}';results[name]={'curves':curves,'cells':selected,'weights':weights.tolist(),'loa_rows':sum(c['loa'] for c in selected),'loa_eligible':int(np.sum(lm>=10.5)) if high else len(lm)}
   ax=axs[row,col]
   for arm,label,color,style in zip(arms+['loa'],['TNG: no brightness cut','Previous colour approximation','DECam integrated r cut','BGS colour/fiber cuts: n=1 proxy','BGS colour/fiber cuts: n=4 proxy','Loa: P12-A / CG15'],['#333333','#888888','#007f98','#c57613','#5a7833','#a42289'],['o--','.:','s-','^--','v:','D-']):
    v=100*np.array(curves[arm]['mean']);ci=100*np.array(curves[arm]['ci95']);ax.plot(range(1,5),v,style,label=label,color=color,lw=1.6)
    if arm in ['bgs_n1','loa']:ax.fill_between(range(1,5),ci[0],ci[1],color=color,alpha=.12)
   ax.set_title(('Shared support' if not high else 'High mass: log M ≥ 10.5')+f' | {width:g}-dex bins\n{results[name]["loa_rows"]:,} Loa galaxies',fontsize=11)
   ax.set_xticks([1,2,3,4],['Lowest','2','3','Highest']);ax.set_ylim(0,100);ax.grid(alpha=.15)
   if col==0:ax.set_ylabel('Quenched fraction (%)')
   if row==1:ax.set_xlabel('Environment quarter within mass–redshift cells')
 fig.suptitle('Improved filter treatment and BGS selection proxies: TNG z=0 versus Loa',fontsize=14)
 fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.955),ncol=3,fontsize=9,frameon=False)
 fig.text(.5,.018,'Same Loa mass–redshift weights across every curve within each panel; Δz = 0.01; shared cells ≥80 galaxies/arm.\nDECam filters integrate a smooth continuum fitted to dust-free TNG griz; no measured spectra or dust model.\nFiber flux uses stellar half-mass radius as a light-radius proxy, Sérsic n=1 or 4, and 1″ seeing.\nR7 truth versus posterior-mean trace ranks; SFR definitions and epoch remain unmatched. Bands: conditional 95% spatial intervals.',ha='center',fontsize=9)
 fig.tight_layout(rect=[0,.12,1,.89]);fig.savefig(O/'refined_overlay.png',dpi=180);fig.savefig(O/'refined_overlay.pdf');plt.close(fig)
 receipt={'results':results,'selection_counts':counts,'checks':checks,'photometry_validation':pval,'source_hashes':{str(p):sha(p) for p in [Path(__file__),O/'inputs.npz',O/'observed_photometry.npz',T/'truth.npz',T/'catalogue.npz',R/'loa_tng_abacus_overlay_20261015/loa_predictions.npz']},'limitations':'No dust, epoch evolution, measured SED or light profile, blend flux, imaging masks/quality, Gaia star separation, fiber assignment or spectroscopic success simulation. Actual desitarget photometric cuts only; all rows assumed real galaxies. Smooth griz continuum extrapolates below g for observed g at z>.2. Fiber sizes are proxies. No P12A validation.'}
 (O/'RESULTS.json').write_text(json.dumps(receipt,indent=2,allow_nan=False))
 print(json.dumps({k:{'rows':v['loa_rows'],'eligible':v['loa_eligible'],'gaps':{a:round(c.get('mean_gap_pp',0),3) for a,c in v['curves'].items()}} for k,v in results.items()},indent=2),flush=True)

if __name__=='__main__':main()
