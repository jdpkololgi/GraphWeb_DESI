"""One-phase pinned altMTL execution; --smoke performs one tracker action."""
import argparse,json,hashlib,subprocess
from pathlib import Path
import multiprocessing
multiprocessing.set_start_method('fork',force=True)
from astropy.table import Table
from LSS.SV3 import altmtltools
import LSS

def main():
 p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--smoke',action='store_true');a=p.parse_args()
 base=a.root.resolve();assert str(base).startswith('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/v1_forward_ph000_20260928/')
 repo=Path(LSS.__file__).resolve().parents[2]
 revision=subprocess.check_output(['git','-c','safe.directory='+str(repo),'-C',str(repo),'rev-parse','HEAD'],text=True).strip()
 if revision!='d942b990860e016e7558c740966fea3e5ce7e7f3':raise RuntimeError('LSS runtime differs from pinned source')
 if not (base/'PREPARATION_READY.json').exists():raise RuntimeError('Preparation validation missing')
 tracker=base/'altmtl0/Univ000/mainsurvey-BRIGHTobscon-TileTracker.ecsv'
 t=Table.read(tracker);before=int((~t['DONEFLAG']).sum());assert before>0
 # Bounded smoke does not need the full target table for an FA-only first action.
 pending=t[~t['DONEFLAG']];pending.sort('ACTIONTIME')
 targets=None if a.smoke and pending[0]['ACTIONTYPE']=='fa' else Table.read(base/'forFA0.fits')
 altmtltools.loop_alt_ledger('BRIGHT',survey='main',zcatdir='/global/cfs/cdirs/desi/spectro/redux/daily/',mtldir='/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/',altmtlbasedir=str(base/'altmtl0'),ndirs=1,debugOrig=True,multiproc=False,nproc=0,targets=targets,mock=True,verbose=True,single_action=a.smoke)
 after=int((~Table.read(tracker)['DONEFLAG']).sum())
 assert after<before,'tracker did not advance'
 (base/('SMOKE.json' if a.smoke else 'ALTMTL_COMPLETE.json')).write_text(json.dumps({'pending_before':before,'pending_after':after,'complete':after==0,'LSS':LSS.__file__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')
 if not a.smoke:assert after==0
if __name__=='__main__':main()
