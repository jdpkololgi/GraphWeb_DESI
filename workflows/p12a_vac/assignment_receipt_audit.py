"""Read-only bounded sample of archived ph006 assignment receipts."""
from pathlib import Path
import json,hashlib,re
import fitsio
BASE=Path('/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl6')
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_photometry_recipe_20260926'
def main():
 dates=sorted((BASE/'Univ000/fa/MAIN').iterdir());chosen=[dates[i] for i in sorted(set([0,len(dates)//4,len(dates)//2,3*len(dates)//4,len(dates)-1]))];records=[]
 for folder in chosen:
  scripts=sorted(folder.glob('fa-*.sh'))
  if not scripts:continue
  p=scripts[0];s=p.read_text();tile=p.stem[3:];f=p.with_name('fba-'+tile+'.fits');h=fitsio.read_header(f,ext=0)
  records.append(dict(script=str(p),output=str(f),script_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),script_text=s,script_mtime=p.stat().st_mtime,output_mtime=f.stat().st_mtime,fa_ver=h.get('FA_VER'),dependencies={k:h[k] for k in h.keys() if k.startswith(('DEPNAM','DEPVER'))},has_fail_fast=('set -e' in s),module_request=re.findall(r'module swap ([^\n]+)',s)))
 src=Path('/global/common/software/desi/perlmutter/desiconda/current/code/LSS/main/py/LSS/main/mockaltmtltools.py');fat=src.parents[1]/'SV3/fatools.py'
 out=dict(records=records,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [src,fat,Path(__file__)]},scope='Five date-spaced first-tile examples, not full FA census. Current LSS code is explanatory evidence, not pinned historical execution. No stdout/stderr receipts found by bounded suffix inventory of altmtl6.',interpretation='Script requests are not verified effective versions. Failed module swap can leave an existing executable and proceed without set -e. Stale scripts, overrides, and reruns also remain possible; no historical failure receipt identifies cause.')
 (OUT/'ASSIGNMENT_RECEIPTS.json').write_text(json.dumps(out,indent=2)+'\n');print([(x['script'],x['module_request'],x['fa_ver']) for x in records])
if __name__=='__main__':main()
