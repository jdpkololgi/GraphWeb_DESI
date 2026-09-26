"""Read back a materialized test catalogue independently of its writer."""
import argparse,json
from pathlib import Path
import fitsio
import numpy as np
from mock_selection_test import save,digest
from verify_mock_test_baseline import OLD,COLS

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);r=p.parse_args().root
 a=json.loads((r/'Y5_MATERIALIZED.json').read_text());idx=np.linspace(0,a['bright']-1,10000,dtype='i8')
 x=fitsio.read(a['path'],rows=idx,columns=COLS);y=fitsio.read(OLD,rows=idx,columns=COLS)
 checks={k:bool(np.array_equal(x[k],y[k])) for k in COLS}
 if not all(checks.values()):raise ValueError('sampled bright content differs')
 with fitsio.FITS(a['path']) as f:
  if f[1].get_nrows()!=a['rows']:raise ValueError('row count changed')
  last=int(f[1].read(rows=[a['rows']-1],columns=['TARGETID'])['TARGETID'][0])
  if last!=a['rows']:raise ValueError('last TARGETID mismatch')
 save(r/'MATERIALIZED_VALIDATION.json',dict(sampled_bright_rows=len(idx),columns_exact=checks,total_rows=a['rows'],last_targetid=last,test_only=True,script_sha256=digest(__file__)))
 print('Readback passes',len(idx),a['rows'],flush=True)
if __name__=='__main__':main()
