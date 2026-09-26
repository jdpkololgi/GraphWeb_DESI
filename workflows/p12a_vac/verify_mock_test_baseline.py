"""Compare regenerated and official ph006 BRIGHT forFA populations, not FAINT RNG."""
import argparse
import json
import os
from pathlib import Path
import fitsio
import numpy as np
from mock_selection_test import digest, save, source

OLD=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/forFA6_nomask.fits')
COLS=['TARGETID','RA','DEC','RSDZ','TRUEZ','R_MAG_APP','R_MAG_ABS','G_R_OBS','BGS_TARGET']

def read(path):
    chunks=[];total=0
    with fitsio.FITS(path) as f:
        for start in range(0,f[1].get_nrows(),500000):
            d=f[1].read(rows=np.arange(start,min(start+500000,f[1].get_nrows())),columns=COLS)
            total+=len(d);chunks.append(d[(d['BGS_TARGET']&2)!=0])
    return np.concatenate(chunks),total

def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);args=p.parse_args();root=args.root
    new=root/'upstream_baseline/SecondGenMocks/AbacusSummitBGS_v2/forFA6_nomask.fits'
    a,na=read(OLD);b,nb=read(new);checks={k:bool(np.array_equal(a[k],b[k])) for k in COLS}
    out=dict(old=source(OLD),new=source(new),old_total=na,new_total=nb,old_bright=len(a),new_bright=len(b),
             exact_ordered_columns=checks,all_bright_columns_exact=all(checks.values()),seed=62026,
             script_sha256=digest(__file__),generated_sha256=digest(new),job=os.environ.get('SLURM_JOB_ID'),
             note='Checks BRIGHT population/IDs exactly; random FAINT downsampling, HIP and subpriority need not reproduce absent historical RNG state. No new fibre assignment run.')
    save(root/'BASELINE_REPLAY.json',out);print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
