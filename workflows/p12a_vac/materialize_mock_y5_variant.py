"""Materialize the verified Y5 test baseline, preserving seeded FAINT realization."""
import argparse,json,os
from pathlib import Path
import fitsio
import numpy as np
from mock_selection_test import digest,save

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);r=p.parse_args().root
 replay=json.loads((r/'Y5_REPLAY.json').read_text())
 if not replay['exact_after_y5_and_reindex']:raise RuntimeError('Y5 replay did not pass')
 baseline=json.loads((r/'BASELINE_REPLAY.json').read_text());src=Path(baseline['new']['path'])
 if digest(src)!=baseline['generated_sha256']:raise RuntimeError('baseline changed')
 out=r/'y5_baseline/forFA6_nomask.fits';out.parent.mkdir(exist_ok=True)
 if out.exists():raise FileExistsError(out)
 count=bright=0
 with fitsio.FITS(src) as f:
  for start in range(0,f[1].get_nrows(),500000):
   d=f[1][start:min(start+500000,f[1].get_nrows())];d=d[d['IN_Y5']!=0]
   d['TARGETID']=np.arange(count+1,count+len(d)+1);count+=len(d);bright+=int(((d['BGS_TARGET']&2)!=0).sum())
   if not out.exists():fitsio.write(out,d,extname='TARGETS',header={'TESTONLY':True,'Y5CUT':True,'SEED':62026,'OBSCON':'BRIGHT'})
   else:
    with fitsio.FITS(out,'rw') as w:w[1].append(d)
 if bright!=replay['old_bright']:raise RuntimeError('unexpected bright count')
 save(r/'Y5_MATERIALIZED.json',dict(path=str(out),rows=count,bright=bright,sha256=digest(out),source_sha256=baseline['generated_sha256'],script_sha256=digest(__file__),job=os.environ['SLURM_JOB_ID'],
    selection='IN_Y5 != 0 applied to pinned seeded upstream preparation; sequential TARGETIDs reassigned after cuts',
    limitations='BRIGHT replay exact; FAINT uses new seed and is not claimed identical to historical stochastic selection. No new fibre assignment or science qualification.'))
 print('MATERIALIZED',count,bright,flush=True)
if __name__=='__main__':main()
