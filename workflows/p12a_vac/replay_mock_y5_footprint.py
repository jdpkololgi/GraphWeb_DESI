"""Test the explicit Y5 footprint toggle against archived BRIGHT preparation."""
import argparse,os
from pathlib import Path
import numpy as np
import healpy as hp
import verify_mock_test_baseline as v
from mock_selection_test import save,digest,source,counts,cap

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);r=p.parse_args().root
 v.COLS.append('IN_Y5')
 old,_=v.read(v.OLD);new,_=v.read(r/'upstream_baseline/SecondGenMocks/AbacusSummitBGS_v2/forFA6_nomask.fits')
 keep=new['IN_Y5']!=0;trimmed=new[keep]
 checks={k:bool(np.array_equal(old[k],trimmed[k])) for k in v.COLS if k!='TARGETID'}
 # The official script assigns sequential IDs after target and footprint cuts.
 reindexed=np.arange(1,len(trimmed)+1,dtype='i8')
 checks['TARGETID_after_postcut_reindex']=bool(np.array_equal(old['TARGETID'],reindexed))
 common=np.load(r/'common.npy');d=new[~keep];k=common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
 extra=counts(d['RSDZ'],cap(d['RA'],d['DEC']),k)
 out=dict(old_bright=len(old),new_bright=len(new),new_bright_inside_y5=int(keep.sum()),old_bright_outside_y5=int((old['IN_Y5']==0).sum()),new_bright_outside_y5=int((~keep).sum()),checks=checks,
          exact_after_y5_and_reindex=all(checks.values()),extra_common_sky_counts=extra.tolist(),script_sha256=digest(__file__),job=os.environ['SLURM_JOB_ID'],
          note='An explicit Y5 cut plus documented postcut ID assignment replays the archived BRIGHT population if all checks pass. This does not identify the historical command or reproduce FAINT RNG/fibre assignment.')
 save(r/'Y5_REPLAY.json',out);print(out,flush=True)
if __name__=='__main__':main()
