"""Recover true mock z for assigned and unassigned targets, by TARGETID."""
from pathlib import Path
import os,json
import numpy as np,fitsio,h5py,hdf5plugin,healpy as hp
import loa_trial as t
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_closure_20260926'
BASE=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2')
def main():
 root=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1');common=np.load(root/'galaxy_angular_support.npy')
 for i in range(2,7):common &=np.load(t.B/f'ph{i:03d}/p3_fields/angular_support_nside256.npz')['support']
 chunks=[];stage={};edges=np.linspace(.15,.55,41)
 def count(name,z,ra,dec,keep):
  cap=t.galactic_cap(ra,dec);keep=keep&common[hp.ang2pix(256,ra,dec,lonlat=True)];stage.setdefault(name,np.zeros((2,40),dtype='i8'))
  for c in [0,1]:stage[name][c]+=np.histogram(z[keep&(cap==c)],edges)[0]
 with fitsio.FITS(BASE/'forFA6_nomask.fits') as f:
  for start in range(0,f[1].get_nrows(),250000):
   d=f[1][start:min(start+250000,f[1].get_nrows())];d=d[(d['BGS_TARGET']&2)!=0];assert np.all(d['R_MAG_APP']<19.5);chunks.append(np.array(d[['TARGETID','RSDZ']]));count('forFA_bright',d['RSDZ'],d['RA'],d['DEC'],np.ones(len(d),bool))
 p=np.concatenate(chunks);p.sort(order='TARGETID');assert np.all(np.diff(p['TARGETID'])>0);del chunks
 for key,path in [('kibo',BASE/'altmtl6/kibo-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.fits'),('loa',BASE/'altmtl6/loa-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.h5')]:
  h5=path.suffix=='.h5';f=h5py.File(path) if h5 else fitsio.FITS(path);tab=f['LSS'] if h5 else f[1];n=len(tab['TARGETID']) if h5 else tab.get_nrows()
  for start in range(0,n,250000):
   sl=slice(start,min(start+250000,n));d={k:tab[k][sl] for k in ['TARGETID','RA','DEC','ZWARN','LOCATION_ASSIGNED','GOODHARDLOC','GOODPRI']} if h5 else tab[sl];at=np.searchsorted(p['TARGETID'],d['TARGETID']);assert np.all(at<len(p));assert np.array_equal(p['TARGETID'][at],d['TARGETID']);z=p['RSDZ'][at]
   for label,keep in [('all',np.ones(len(z),bool)),('assigned',d['LOCATION_ASSIGNED'].astype(bool)),('zwarn0',d['ZWARN']==0),('usable',d['GOODHARDLOC'].astype(bool)&d['GOODPRI'].astype(bool))]:count(key+'_'+label,z,d['RA'],d['DEC'],keep)
  f.close()
 out=dict(stage={k:v.tolist() for k,v in stage.items()},forfa_bright_rows=len(p),all_ids_joined=True,all_bright_rlt19p5=True,job=os.environ.get('SLURM_JOB_ID'),script_sha256=t.digest(__file__),note='True RSDZ from forFA for all targets, including unassigned; fixed common sky. No repairs. Stage survival is for the mock population, not a measurement of DESI true-z completeness.')
 t.save(OUT/'MOCK_RETENTION.json',out)
 for k,v in stage.items():print(k,v.reshape(2,4,10).sum(-1).tolist(),flush=True)
if __name__=='__main__':main()
