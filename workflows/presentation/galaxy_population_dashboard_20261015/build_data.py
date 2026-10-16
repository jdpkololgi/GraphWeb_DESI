"""Build aggregate-only Page data; no object identifiers leave the local cache."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import sys,json,hashlib,gzip,base64,struct
import numpy as np
import fitsio
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
ROOT=Path('/global/u2/d/dkololgi');R=ROOT/'fmsu_tokyo_20261010';O=R/'galaxy_population_dashboard_20261015';O.mkdir(exist_ok=True)
VAC=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')
PHOT=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full_HPmapcut.dat.fits')
MB=np.arange(9,12.0001,.5);RB=np.array([-100,16,17,18,18.5,19,19.5,20,20.5,100])
ZB=np.array([.01,.02,.03,.04,.05,.06,.15,.20,.25,.30,.35,.40,.45,.50,.55])
TZ=np.array([.025,.05,.175,.225,.25,.275,.325,.375,.425,.475,.525])
PROPS=['quenched','mass','sfr','ssfr','colour','age','av','fagn']
checks=[]

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

def obs_props(a):
 m=np.log10(np.maximum(a['MASS_CG_15'],1e-30));s=np.log10(np.maximum(a['SFR_CG_15'],1e-30));age=np.log10(np.maximum(a['AGE_CG_15'],1e-30)/1000)
 good=(m>=9)&(m<12)
 for k in ['MASS_CG_15','SFR_CG_15','AGE_CG_15']:good&=np.isfinite(a[k])&(a[k]>0)
 good&=np.isfinite(a['AV_CG_15'])&(a['AV_CG_15']>=0)&np.isfinite(a['FAGN_CG_15'])&(a['FAGN_CG_15']>=0)&(a['FAGN_CG_15']<=1)
 y=np.column_stack([s-m<-11,m,np.maximum(s,-14),np.maximum(s-m,-14),np.zeros(len(m)),age,a['AV_CG_15'],a['FAGN_CG_15']])
 return good,y

def current():
 cache=O/'current.npz'
 if cache.exists():return dict(np.load(cache))
 cols=['TARGETID','RA','DEC','Z','SUPPORTED','CIGALE_MATCHED','CIGALE_AMBIGUOUS','CIGALE_SPECTYPE','MASS_CG_15','SFR_CG_15','AGE_CG_15','AV_CG_15','FAGN_CG_15','P_WEB','EIGENVALUE_MEAN']
 out=[]
 with fitsio.FITS(VAC) as f:
  n=f[1].get_nrows()
  for start in range(0,n,200000):
   a=f[1].read(rows=np.arange(start,min(start+200000,n)),columns=cols)
   good,y=obs_props(a);good&=a['SUPPORTED'].astype(bool)&a['CIGALE_MATCHED'].astype(bool)&~a['CIGALE_AMBIGUOUS'].astype(bool)&(np.char.strip(a['CIGALE_SPECTYPE'].astype('U'))=='GALAXY')&(a['Z']>=.15)&(a['Z']<.55)
   out.append(a[good])
   if start%1000000==0:print('current scan',start,n,flush=True)
 a=np.concatenate(out);del out;good,y=obs_props(a)
 order=np.argsort(a['TARGETID']);sid=a['TARGETID'][order];assert np.all(np.diff(sid)>0)
 r=np.full(len(a),np.nan);color=r.copy();found=np.zeros(len(a),np.uint8)
 with fitsio.FITS(PHOT) as f:
  n=f[1].get_nrows()
  for start in range(0,n,200000):
   b=f[1].read(rows=np.arange(start,min(start+200000,n)),columns=['TARGETID','RA','DEC','Z_not4clus','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R','PHOTSYS'])
   at=np.searchsorted(sid,b['TARGETID']);ok=(at<len(sid))&(sid[np.minimum(at,len(sid)-1)]==b['TARGETID']);ii=order[at[ok]];b=b[ok]
   for x,k in [('RA','RA'),('DEC','DEC'),('Z_not4clus','Z')]:assert np.array_equal(b[x],a[k][ii])
   fr=b['FLUX_R']/b['MW_TRANSMISSION_R'];fg=b['FLUX_G']/b['MW_TRANSMISSION_G'];ps=np.char.strip(b['PHOTSYS'].astype('U'));assert np.isin(ps,['N','S']).all()
   valid=np.isfinite(fr)&(fr>0)&np.isfinite(fg)&(fg>0)
   r[ii[valid]]=22.5-2.5*np.log10(fr[valid])-np.where(ps[valid]=='N',.04,0)
   color[ii[valid]]=-2.5*np.log10(fg[valid]/fr[valid]);found[ii]+=1
   if start%2000000==0:print('photometry scan',start,n,flush=True)
 assert np.all(found==1);valid=np.isfinite(r)&np.isfinite(color)&(r<19.5)&np.isfinite(a['EIGENVALUE_MEAN']).all(1)
 y[:,4]=color
 d={'y':y[valid].astype('f4'),'z':a['Z'][valid].astype('f4'),'r':r[valid].astype('f4'),'p':a['P_WEB'][valid].astype('f4'),'trace':a['EIGENVALUE_MEAN'][valid].sum(1).astype('f4')}
 assert np.allclose(d['p'].sum(1),1,atol=1e-5)
 np.savez_compressed(cache,**d);print('current cache',len(d['z']),flush=True);return d

def legacy():
 a=dict(np.load(R/'legacy_gat_comparison_20261014/legacy_cigale.npz'));good,y=obs_props(a)
 good&=(np.char.strip(a['SPECTYPE'].astype('U'))=='GALAXY')&(a['z']>=.01)&(a['z']<.06)&(a['FLUX_G']>0)&(a['FLUX_R']>0)
 # Preserve the saved CG15-associated flux convention, do not silently call it LS dust-corrected.
 r=22.5-2.5*np.log10(np.maximum(a['FLUX_R'],1e-30));y[:,4]=-2.5*np.log10(np.maximum(a['FLUX_G'],1e-30)/np.maximum(a['FLUX_R'],1e-30))
 return {'y':y[good],'z':a['z'][good],'r':r[good],'p':a['P_WEB'][good]}

def tng():
 t=R/'tng_density_control_20261015';d=np.load(t/'catalogue.npz');p=np.load(R/'tng_bright_refined_20261015/inputs.npz');ids=p['ids'];m=np.log10(d['mass_msun'][ids]);s=np.log10(np.maximum(d['sfr'][ids],1e-14));ss=np.log10(np.maximum(d['sfr'][ids]/d['mass_msun'][ids],1e-14))
 y=np.column_stack([d['sfr'][ids]/d['mass_msun'][ids]<1e-11,m,s,ss,d['g_r'][ids],np.zeros((len(ids),3))])
 cache=O/'tng_observing.npz'
 if cache.exists():obs=np.load(cache)['gr']
 else:
  sys.path.insert(0,'/global/common/software/desi/perlmutter/desiconda/current/code/speclite/v1.0.0/lib/python3.13/site-packages')
  from speclite import filters
  wave=np.arange(2900.,11501.,10.);fs=filters.load_filters('sdss2010-g','sdss2010-r','sdss2010-i','sdss2010-z');knots=np.array([f.effective_wavelength.value for f in fs]);basis=interp1d(np.log(knots),np.eye(4),axis=0,fill_value='extrapolate')
  w=[]
  for f in filters.load_filters('decam2014-g','decam2014-r'):
   a=np.interp(wave,f.wavelength,f.response,left=0,right=0)*np.gradient(wave)/wave;w.append(a/a.sum())
  w=np.array(w);nodes=np.load(R/'tng_bright_refined_20261015/observed_photometry.npz')['logfnu_nodes'];obs=np.empty((len(TZ),len(ids),2),dtype='f4')
  for j,z in enumerate(TZ):
   b=basis(np.log(wave/(1+z))).T;dm=Planck18.distmod(z).value
   for st in range(0,len(ids),2048):
    ix=slice(st,st+2048);obs[j,ix]=-2.5*np.log10((10**(nodes[ix]@b)@w.T)*(1+z))+dm
   print('TNG observing redshift',z,flush=True)
  np.savez_compressed(cache,gr=obs)
 eig=np.load(t/'truth.npz');return y,{scale:eig[f'eig_R{scale}_n512'][ids] for scale in [2,7]},obs

# Fixed, parent-defined environment quarters within each mass/redshift cell.
# These stay fixed as users change brightness; no claim of reranked selected samples.
def fractional_quarters(x):
 _,inv,n=np.unique(x,return_inverse=True,return_counts=True);end=np.cumsum(n);start=end-n;N=len(x)
 return np.array([np.maximum(0,np.minimum(end,(j+1)*N/4)-np.maximum(start,j*N/4))/n for j in range(4)]).T[inv]

records=[];audits=[]
def aggregate(source,zi,scale,y,r,p,trace=None):
 mb=np.searchsorted(MB,y[:,1],side='right')-1;rb=np.searchsorted(RB,r,side='right')-1
 assert ((mb>=0)&(mb<6)&(rb>=0)&(rb<9)).all()
 weights={0:np.eye(4)[p.argmax(1)]}
 if source!=2:weights[1]=p
 if trace is not None:
  q=np.zeros_like(p,dtype=float)
  for m in np.unique(mb):
   ii=np.flatnonzero(mb==m);q[ii]=fractional_quarters(trace[ii])
  weights[2]=q
 for mode,w in weights.items():
  cell=mb*9+rb;total=0
  for env in range(4):
   ww=w[:,env];n=np.bincount(cell,weights=ww,minlength=54);n2=np.bincount(cell,weights=ww**2,minlength=54)
   sums=np.array([np.bincount(cell,weights=ww*y[:,j],minlength=54) for j in range(8)]).T
   for c in np.flatnonzero(n>1e-8):
    mean=sums[c]/n[c];assert np.all(np.isfinite(mean));assert np.all(np.abs(mean)<32)
    records.append([source,zi,int(c//9),int(c%9),mode,scale,env,float(n[c]),float(n2[c]),*np.rint(mean*1000).astype(int).tolist()]);total+=n[c]
  assert abs(total-len(y))<max(1e-3,len(y)*1e-5)
  audits.append({'source':source,'z_index':zi,'scale':scale,'method':mode,'rows':len(y),'aggregate_weight':total})

def main():
 for source,get in [(0,legacy),(1,current)]:
  d=get();zz=np.searchsorted(ZB,d['z'],side='right')-1
  for zi in np.unique(zz):
   assert zi!=5 # no interpolation or rows in the VAC gap
   ii=zz==zi;aggregate(source,int(zi),0,d['y'][ii],d['r'][ii],d['p'][ii],d.get('trace',None)[ii] if 'trace' in d else None)
  print('aggregated source',source,len(d['z']),len(records),flush=True)
 y,eig,obs=tng()
 for zi,z in enumerate(TZ):
  yy=y.copy();yy[:,4]=obs[zi,:,0]-obs[zi,:,1]
  for scale,e in eig.items():
   p=np.eye(4)[(e>.2).sum(1)];aggregate(2,zi,scale,yy,obs[zi,:,1],p,e.sum(1))
 print('total records',len(records),flush=True)
 # Little-endian fixed records: 7 category bytes, weight and squared weight,
 # eight property means quantized at .001 native units. No target IDs/coordinates.
 raw=b''.join(struct.pack('<7B2f8h',*row) for row in records)
 packed=base64.b64encode(gzip.compress(raw,compresslevel=9,mtime=0)).decode()
 meta={'properties':PROPS,'massEdges':MB.tolist(),'rEdges':RB.tolist(),'zEdges':ZB.tolist(),'tngObservingZ':TZ.tolist(),'recordBytes':31,'records':len(records),'scale':1000,'data':packed,
 'sources':[{'id':0,'name':'Loa low-z GAT','coverage':[.01,.06],'environment':'Saved GAT native class or class probability; no tidal trace','properties':'CG15; saved associated fluxes','colour':'Observed g-r from saved CG15 fluxes; flux convention differs from current LS photometry'},
 {'id':1,'name':'Loa P12-A','coverage':[.15,.55],'environment':'P_WEB native posterior class probabilities; trace of posterior-mean R7 eigenvalues','properties':'CG15; supported unique matches; explicit Bright selection','colour':'Milky-Way-corrected observed Legacy Survey g-r; r control uses north r-0.04'},
 {'id':2,'name':'TNG300 z=0','coverage':[0,0],'environment':'R2/R7 simulation truth, eigenvalue threshold0.2','properties':'Held-out galaxies, logM9-11.5; instantaneous gas SFR; other CG15 properties unavailable','colour':'Synthetic observer-frame dust-free DECam g-r; observing redshift is not evolution'}],
 'notes':['All Loa panels share valid CG15 mass/SFR/age/AV/FAGN and photometry eligibility; property changes are means, not medians.','Rank quarters are fixed in the parent mass-redshift bins before brightness selection; legacy GAT trace ranks unavailable.','Native class identities are source-specific; equal labels or ranks do not establish matched physical environments.','Means are quantized to .001 native units; counts and squared weights retain float32 precision.','No inferential uncertainty is claimed for interactive mixtures.','TNG log sSFR and log SFR have -14 floors for zero/very small SFR.','No interpolation across0.06-0.15; any combined sample remains an explicitly disjoint mixture.'],
 'sourceStats':{str(p):{'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in [VAC,PHOT,R/'legacy_gat_comparison_20261014/legacy_cigale.npz',R/'tng_density_control_20261015/catalogue.npz',R/'tng_density_control_20261015/truth.npz']}}
 (O/'data.json').write_text(json.dumps(meta,separators=(',',':'),allow_nan=False));(O/'aggregation_audit.json').write_text(json.dumps(audits,indent=2));print('packed bytes',len(packed),'json', (O/'data.json').stat().st_size,flush=True)
 # Exact local baseline for frontend numerical parity.
 (O/'reference_records.json').write_text(json.dumps(records[:20]));(O/'BUILD.json').write_text(json.dumps({'builder_sha256':sha(__file__),'records':len(records),'uncompressed_bytes':len(raw),'data_sha256':sha(O/'data.json'),'audited_aggregations':len(audits)},indent=2))

if __name__=='__main__':main()
