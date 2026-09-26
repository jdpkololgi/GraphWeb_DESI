"""Separate changed footprint/rows from changed photometry in LSS baseline replay."""
import argparse,json,os
import numpy as np
import healpy as hp
from verify_mock_test_baseline import read,OLD,COLS
from mock_selection_test import save,source,digest,cap,counts

def main():
 from pathlib import Path
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);r=p.parse_args().root
 new=r/'upstream_baseline/SecondGenMocks/AbacusSummitBGS_v2/forFA6_nomask.fits'
 a,_=read(OLD);b,_=read(new)
 # Use sky coordinates plus redshift; fail rather than collapse repeated directions.
 def keys(d):
  k=np.empty(len(d),dtype=[('RA','f8'),('DEC','f8'),('RSDZ','f8')]);k['RA']=d['RA'];k['DEC']=d['DEC'];k['RSDZ']=d['RSDZ'];return k
 ka=keys(a);kb=keys(b);order=np.argsort(ka,order=['RA','DEC','RSDZ']);sa=ka[order]
 if np.any(sa[1:]==sa[:-1]):raise ValueError('duplicate old 3D positions require a richer join')
 bo=np.argsort(kb,order=['RA','DEC','RSDZ']);sb=kb[bo]
 if np.any(sb[1:]==sb[:-1]):raise ValueError('duplicate new 3D positions require a richer join')
 at=np.searchsorted(sa,kb);valid=at<len(sa);matched=np.zeros(len(b),bool);matched[valid]=sa[at[valid]]==kb[valid]
 ia=order[at[matched]];ib=np.flatnonzero(matched);oldonly=np.ones(len(a),bool);oldonly[ia]=False
 checks={k:dict(exact=bool(np.array_equal(a[k][ia],b[k][ib])),max_abs_difference=float(np.max(np.abs(a[k][ia]-b[k][ib])))) for k in COLS if k not in ['RA','DEC']}
 common=np.load(r/'common.npy');stage={}
 for name,d in [('old',a),('new',b),('old_only',a[oldonly]),('new_only',b[~matched])]:
  pix=hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True);stage[name]=counts(d['RSDZ'],cap(d['RA'],d['DEC']),common[pix]).tolist()
 out=dict(old=source(OLD),new=source(new),shared=len(ia),old_only=int(oldonly.sum()),new_only=int((~matched).sum()),shared_columns=checks,common_sky_counts=stage,
          script_sha256=digest(__file__),job=os.environ['SLURM_JOB_ID'],note='Exact RA/DEC/RSDZ join; changed TARGETID expected when footprint changes ordering. Differing footprint requires historical tile list to attribute, not inferred from filename.')
 save(r/'BASELINE_JOIN.json',out);print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
