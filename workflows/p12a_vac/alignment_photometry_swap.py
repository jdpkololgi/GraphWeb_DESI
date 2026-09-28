"""Exploratory fixed-row mapping sensitivity; not a physical generator replacement."""
import sys,json,os,hashlib
from pathlib import Path
import numpy as np
import fitsio,healpy as hp
from alignment_parent_screen import BASE,MASK
from mock_selection_test import cap
from photometry_recipe_fingerprint import k
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'docs/evidence/p12a_canonical_v1_replay_20260928/source'))
from hodpy.k_correction import DESI_KCorrection

def main():
 out={'job':os.environ.get('SLURM_JOB_ID'),'stride':100,'warning':'conditional mapping swap on delivered flux-limited parents; conventions and parent truncation preclude causal HOD attribution','versions':{}}
 common=np.load(MASK); edges=np.array([.15,.25,.35,.45,.55]); ks={p:DESI_KCorrection('r',p) for p in ['N','S']}
 for v in ['v0.1','v1']:
  path=BASE/v/'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
  with fitsio.FITS(path) as f:
   n=f[1].get_nrows(); d=f[1].read(rows=np.arange(0,n,100),columns=['RA','DEC','Z','R_MAG_APP','G_R_REST'])
  z=d['Z'];d=d[(z>=.15)&(z<.55)];d=d[common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]]
  z=d['Z'].astype(float);c=d['G_R_REST'].astype(float);north=d['DEC']>32.375
  new=np.where(north,ks['N'].k(z,c),ks['S'].k(z,c))-.67*(z-.1)
  old=k(z,c,'r')-.8*(z-.1);delta=new-old
  ca=cap(d['RA'],d['DEC']);r=d['R_MAG_APP']
  variants={}
  for name,mag in [('stored',r),('swapped',r+delta if v=='v0.1' else r-delta)]:
   variants[name]=[np.histogram(z[(ca==i)&(mag>=12)&(mag<19.5)],edges)[0].tolist() for i in [0,1]]
  out['versions'][v]={'path':str(path),'n_total':n,'n_sample_common':len(d),'counts_sample':variants,'delta_new_minus_old_median_by_z':[float(np.median(delta[(z>=lo)&(z<hi)])) for lo,hi in zip(edges[:-1],edges[1:])]}
 loa=np.load(ROOT/'docs/evidence/p12a_joint_flags_20260926/HISTOGRAMS.npz')['loa'][:,:,:75].sum(-1).reshape(2,4,10).sum(-1)
 out['loa_observed_counts']=loa.tolist()
 out['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 target=ROOT/'docs/evidence/p12a_photometry_swap_20260928';target.mkdir(exist_ok=True)
 (target/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
