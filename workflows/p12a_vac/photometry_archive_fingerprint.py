"""Compare mapping identities in three existing ph006 raw revisions; no joins."""
from pathlib import Path
import os,json,hashlib
import numpy as np,fitsio
from astropy.cosmology import FlatLambdaCDM
from photometry_recipe_fingerprint import k,stats,OUT

def main():
 base=Path('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS');out={}
 grid=np.linspace(0,1,20001);cosmo=FlatLambdaCDM(H0=67.36,Om0=(.02237+.12+.0006442)/.6736**2,Tcmb0=0);dl=cosmo.luminosity_distance(grid).value*cosmo.h
 for rev in ['v0/z0.200','v0.1/z0.200/old','v0.1/z0.200']:
  p=base/rev/'cutsky_BGS_z0.200_AbacusSummit_base_c000_ph006.fits'
  with fitsio.FITS(p) as f:
   n=f[1].get_nrows();d=f[1].read(rows=np.linspace(0,n-1,10000,dtype=int),columns=['Z','Z_COSMO','R_MAG_APP','R_MAG_ABS','G_R_REST','G_R_OBS'])
  z=d['Z'];c=d['G_R_REST'];sel=(z>=.15)&(z<.55)&(d['R_MAG_APP']<19.5);kr=k(z,c,'r');dc=d['G_R_OBS']-c-k(z,c,'g')+kr
  dr=d['R_MAG_APP']-d['R_MAG_ABS']-5*np.log10(np.maximum(np.interp(z,grid,dl),1e-20))-25-kr
  fit=np.polyfit(z[sel]-.1,dr[sel],1)
  distance_variants={}
  for zn in ['Z','Z_COSMO']:
   dz=d[zn]; rr=d['R_MAG_APP']-d['R_MAG_ABS']-5*np.log10(np.maximum(np.interp(dz,grid,dl),1e-20))-25-kr;ff=np.polyfit(z[sel]-.1,rr[sel],1)
   distance_variants[zn]=dict(fit=ff.tolist(),after_fit=stats((rr-np.polyval(ff,z-.1))[sel]))
  out[rev]=dict(distance_variants=distance_variants,path=str(p),total_rows=n,sample_rows=len(d),selected=int(sel.sum()),colour_residual=stats(dc[sel]),r_fit_slope_intercept=fit.tolist(),r_after_fit=stats((dr-np.polyval(fit,z-.1))[sel]))
 result=dict(revisions=out,job=os.environ.get('SLURM_JOB_ID'),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),note='Uniformly spaced 10000 raw rows per revision; independent samples, not galaxy joins or matched count comparisons. Approximate native distance as main fingerprint. No version is a replacement candidate.')
 (OUT/'ARCHIVE_FINGERPRINT.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
