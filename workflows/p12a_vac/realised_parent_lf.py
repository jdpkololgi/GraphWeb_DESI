"""Delivered-parent luminosity diagnostic; no assumption of parent completeness."""
import json,os,hashlib
from pathlib import Path
import numpy as np,healpy as hp
from astropy.io import fits
from scipy.integrate import cumulative_trapezoid
from alignment_parent_screen import BASE,MASK
R=Path(__file__).resolve().parents[2];O=R/'docs/evidence/p12a_realised_lf_20260928';O.mkdir(exist_ok=True)
mask=np.load(MASK);ze=np.linspace(.15,.55,5);limits=np.array([-21.,-22.,-23.,-24.]);grid=np.linspace(0,1,100001);om=(.02237+.1200+.00064420)/.6736**2;chi=cumulative_trapezoid(299792.458/100/np.sqrt(om*(1+grid)**3+1-om),grid,initial=0)
volume=mask.sum()*hp.nside2pixarea(256)*np.diff(np.interp(ze,grid,chi)**3)/3
out={}
for version in ['v0.1','v1']:
 p=BASE/version/'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits';stat=[p.stat().st_size,p.stat().st_mtime_ns];h=np.zeros((4,4));bright=h.copy();red=h.copy();rmax=np.full((4,4),-np.inf)
 with fits.open(p,memmap=True) as f:
  for start in range(0,len(f[1].data),2000000):
   d={c:np.array(f[1].data[c][start:start+2000000:20]) for c in ['RA','DEC','Z','R_MAG_ABS','R_MAG_APP','G_R_REST']}
   z=d['Z'];sel=mask[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]&(z>=.15)&(z<.55)
   assert np.all(np.isfinite(d['R_MAG_ABS'][sel]))
   for j,m in enumerate(limits):
    take=sel&(d['R_MAG_ABS']<m);h[:,j]+=np.histogram(z[take],ze)[0]*20
    bright[:,j]+=np.histogram(z[take&(d['R_MAG_APP']>=12)&(d['R_MAG_APP']<19.5)],ze)[0]*20
    red[:,j]+=np.histogram(z[take&(d['G_R_REST']>=.75)],ze)[0]*20
    for i in range(4):
     rr=d['R_MAG_APP'][take&(z>=ze[i])&(z<ze[i+1])]
     if len(rr):rmax[i,j]=max(rmax[i,j],float(rr.max()))
 assert stat==[p.stat().st_size,p.stat().st_mtime_ns]
 out[version]={'path':str(p),'stat':stat,'estimated_counts':h.tolist(),'density':(h/volume[:,None]).tolist(),'bright_fraction':np.divide(bright,h,out=np.zeros_like(h),where=h>0).tolist(),'red_fraction':np.divide(red,h,out=np.zeros_like(h),where=h>0).tolist(),'maximum_delivered_r':np.where(np.isfinite(rmax),rmax,np.nan).tolist()}
lf={}
for b in ['wsys','nowsys']:
 p=R/f'docs/evidence/p12a_photometry_swap_20260928/lf_{b}.dat';x,y=np.loadtxt(p,unpack=True);n=10**np.interp(limits,x,y);lf[b]={'cumulative':n.tolist(),'v1_density_over_table':(np.array(out['v1']['density'])/n).tolist()}
result={'parents':out,'lf':lf,'M_limits':limits.tolist(),'z_edges':ze.tolist(),'volume_Mpch3':volume.tolist(),'job':os.environ.get('SLURM_JOB_ID'),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'limits':'Stored R_MAG_ABS; delivered parent density is a lower bound under missing faint support, NOT a complete LF estimate. Numerical reference convention must be checked. Common geometric volume, approximate c000 distance, RSD shells, deterministic stride20 ph000; no errors/covariance. No mutation.'}
(O/'RESULTS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
