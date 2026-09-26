"""Bounded ph006 raw-cutsky census on the same matched sky as count plots."""
import json
from pathlib import Path
import numpy as np,fitsio,healpy as hp
import loa_trial as t
ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1')
OUT=Path(__file__).resolve().parents[2]/'docs/figures/p12a_selection_20260925'
def main():
 common=np.load(ROOT/'galaxy_angular_support.npy')
 for i in range(2,7):common &=np.load(t.B/f'ph{i:03d}/p3_fields/angular_support_nside256.npz')['support']
 raw=Path('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v0.1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph006.fits');edges=np.linspace(.15,.60,46);hist={};abedges=np.linspace(-27,-17,101);abh=np.zeros((45,100),dtype='i8');ranges={};total=0
 def add(label,z,cap,k):
  hist.setdefault(label,np.zeros((2,45),dtype='i8'))
  for c in [0,1]:hist[label][c]+=np.histogram(z[k&(cap==c)],edges)[0]
 with fitsio.FITS(raw) as f:
  columns=f[1].get_colnames();n=f[1].get_nrows()
  for start in range(0,n,500000):
   d=f[1][start:min(start+500000,n)];total+=len(d);pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);k=common[pix];d=d[k];cap=t.galactic_cap(d['RA'],d['DEC']);z=d['Z'];m=d['R_MAG_APP'];bright=m<19.5
   for limit in [19.4,19.5,19.54,19.6,19.7,20.0]:add(f'raw_rlt{limit}',z,cap,m<limit)
   add('raw_rlt19.5_zcosmo',d['Z_COSMO'],cap,bright);add('raw_rlt19.5_y5',z,cap,bright&(d['IN_Y5']!=0))
   abh+=np.histogram2d(z[bright],d['R_MAG_ABS'][bright],bins=[edges,abedges])[0].astype('i8')
   for col in ['Z','Z_COSMO','R_MAG_APP','R_MAG_ABS']:
    if len(d):
     v=d[col];r=ranges.setdefault(col,[float('inf'),-float('inf')]);r[0]=min(r[0],float(v.min()));r[1]=max(r[1],float(v.max()))
   if start%10000000==0:print('raw',start,n,flush=True)
 out=dict(path=str(raw),bytes=raw.stat().st_size,mtime_ns=raw.stat().st_mtime_ns,rows=total,columns=columns,edges=edges.tolist(),counts={k:v.tolist() for k,v in hist.items()},absolute_magnitude_edges=abedges.tolist(),absolute_magnitude_histogram=abh.tolist(),ranges_common_sky=ranges,note='Raw parent before fibre assignment, same common pixels. Magnitude ablations are diagnostic, not candidate cuts. Z and Z_COSMO separated. No HOD or luminosity-evolution cause inferred from n(z) alone.',script_sha256=t.digest(__file__))
 t.save(OUT/'UPSTREAM.json',out)
 for k,v in hist.items():print(k,v[:,:40].reshape(2,4,10).sum(-1).tolist(),flush=True)
if __name__=='__main__':main()
