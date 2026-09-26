"""Diagnostic paired SDSS/DECam synthesis for existing Loa targets.

Uses pinned FastSpecFit 3.1.5 template implementation and speclite 0.20.
FastSpecFit output FLUX columns are already de-reddened (3.1.5 fastspecfit.py:99-105).
Broadband differences only; does not synthesize a Petrosian aperture or repair data.
"""
from pathlib import Path
import sys,os,json,hashlib
ROOT=Path('/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/code')
for p in [ROOT/'fastspecfit/3.1.5/lib/python3.10/site-packages',ROOT/'speclite/v0.20/lib/python3.10/site-packages',ROOT/'desiutil/3.4.3/lib/python3.10/site-packages']:
 sys.path.insert(0,str(p))
import numpy as np,fitsio,healpy as hp
from fastspecfit.templates import Templates
from speclite.filters import load_filters
from photometry_recipe_fingerprint import k
import loa_trial as t
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_passband_followup_20260926'
CAT=Path('/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/fastspec-loa-main-bright-nside1-hp00.fits')
TPL=Path('/global/cfs/cdirs/desi/public/external/templates/fastspecfit/2.0.0/ftemplates-chabrier-2.0.0.fits')
def quant(a):return dict(n=len(a),median=float(np.median(a)),p16_p84=np.quantile(a,[.16,.84]).tolist(),max_abs=float(np.max(np.abs(a))))
def main():
 root=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
 with np.load(root/'catalogue.npz') as c:ids=np.sort(c['targetid'][(c['shell']>=0)&(c['cap']==0)])
 common=np.load(root/'galaxy_angular_support.npy')
 for i in range(2,7):common &=np.load(t.B/f'ph{i:03d}/p3_fields/angular_support_nside256.npz')['support']
 with fitsio.FITS(CAT) as f:
  m=f['METADATA'].read(columns=['TARGETID','Z','RA','DEC','PHOTSYS','FLUX_G','FLUX_R','MW_TRANSMISSION_G','MW_TRANSMISSION_R'])
  at=np.searchsorted(ids,m['TARGETID']);joined=(at<len(ids))&(ids[np.minimum(at,len(ids)-1)]==m['TARGETID']);r=22.5-2.5*np.log10(np.maximum(m['FLUX_R'],1e-20))
  keep=joined&(m['Z']>=.15)&(m['Z']<.55)&(r>=12)&(r<19.5)&(np.char.strip(m['PHOTSYS'].astype('U'))=='S')&common[hp.ang2pix(256,m['RA'],m['DEC'],lonlat=True)]
  rng=np.random.default_rng(20260926);rows=[];eligible=[]
  for lo,hi in zip(np.linspace(.15,.55,17)[:-1],np.linspace(.15,.55,17)[1:]):
   ix=np.flatnonzero(keep&(m['Z']>=lo)&(m['Z']<hi));eligible.append(len(ix));rows.extend(rng.choice(ix,min(32,len(ix)),replace=False).tolist())
  rows=np.sort(rows);m=m[rows];s=f['SPECPHOT'].read(rows=rows,columns=['TARGETID','COEFF','TAUV','VDISP','RCHI2_PHOT','FLUX_SYNTH_PHOTMODEL_G','FLUX_SYNTH_PHOTMODEL_R','ABSMAG01_SYNTH_SDSS_G','ABSMAG01_SYNTH_SDSS_R'])
  assert np.array_equal(m['TARGETID'],s['TARGETID']);header=str(f[0].read_header())
 print('selected',len(rows),'eligible',eligible,flush=True)
 tpl=Templates(template_file=str(TPL));assert s['COEFF'].shape[1]==tpl.ntemplates
 filters=load_filters('decam2014-g','decam2014-r','sdss2010-g','sdss2010-r');restfilters=[filters[i].create_shifted(.1) for i in [2,3]]
 values=[];valid=[]
 for i,(mm,ss) in enumerate(zip(m,s)):
  if not np.all(np.isfinite(ss['COEFF'])) or np.sum(ss['COEFF'])<=0 or not np.isfinite(ss['TAUV']):continue
  z=float(mm['Z']);flux=tpl.convolve_vdisp(ss['COEFF'].dot(tpl.flux),float(ss['VDISP']))*np.exp(-ss['TAUV']*tpl.dust_klambda)
  mags=np.array([float(f.get_ab_magnitude(flux,tpl.wave*(1+z))) for f in filters]);restmag=np.array([float(f.get_ab_magnitude(flux,tpl.wave)) for f in restfilters]);restgr=restmag[0]-restmag[1]
  kgkr=mags[2:]-restmag+2.5*np.log10(1+z)
  gama=np.array([k(np.array([z]),np.array([restgr]),band)[0] for band in ['g','r']])
  measured=-2.5*np.log10(mm['FLUX_G']/mm['FLUX_R'])
  published=-2.5*np.log10(ss['FLUX_SYNTH_PHOTMODEL_G']/ss['FLUX_SYNTH_PHOTMODEL_R'])
  # All colour differences in magnitudes; paired filters on precisely the same SED.
  values.append([z,mags[3]-mags[1],(mags[2]-mags[3])-(mags[0]-mags[1]),mags[0]-mags[1]-published,restgr-(ss['ABSMAG01_SYNTH_SDSS_G']-ss['ABSMAG01_SYNTH_SDSS_R']),measured-published,kgkr[1]-gama[1],(kgkr[0]-kgkr[1])-(gama[0]-gama[1]),restgr,float(ss['RCHI2_PHOT'])]);valid.append(i)
 a=np.array(values);names=['z','sdss_minus_decam_r','sdss_minus_decam_gr','replay_decam_gr_error','replay_rest_sdss_gr_error','observed_minus_model_decam_gr','sed_minus_gama_Kr','sed_minus_gama_KgKr','rest_sdss_gr','rchi2_phot']
 np.savez_compressed(OUT/'SAME_GALAXY_SAMPLE.npz',values=a,names=np.array(names),targetid=m['TARGETID'][valid],catalogue_rows=np.array(rows)[valid])
 result=dict(catalogue=str(CAT),template=str(TPL),header=header,job=os.environ.get('SLURM_JOB_ID'),eligible_per_dz025=eligible,selected=len(rows),valid=len(a),shells=[],script_sha256=t.digest(__file__),implementation_hashes={str(p):t.digest(p) for p in [Path(sys.modules['fastspecfit.templates'].__file__),Path(sys.modules['speclite.filters'].__file__)]},caveats='One healpix hp00, southern common-sky VAC targets, 32 random rows per dz=.025 bin before fit validity. Not survey representative. Stellar+nebular templates, fitted dust and dispersion; infrared dust emission and IGM omitted in optical. Replay colours validated against catalogue. No aperture/flux-estimator conversion; no observed SDSS Petrosian measurement. K comparison is fitted-SED dependent, not independent truth.')
 for lo,hi in zip([.15,.25,.35,.45],[.25,.35,.45,.55]):
  b=a[(a[:,0]>=lo)&(a[:,0]<hi)];result['shells'].append(dict(z=[lo,hi],stats={n:quant(b[:,j]) for j,n in enumerate(names) if j and len(b)}))
 result['replay_pass']=bool(np.max(np.abs(a[:,3:5]))<.005)
 (OUT/'SAME_GALAXY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
 if not result['replay_pass']:raise RuntimeError('Reconstruction parity failed; do not interpret passband result until resolved')
if __name__=='__main__':main()
