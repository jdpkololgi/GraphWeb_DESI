"""Read-only conditional comparison of existing magnitude term and paired filters.

This is algebra on saved diagnostics, not a fitted offset or catalogue repair.
"""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'docs/evidence/p12a_nz_magnitude_20260926'

def main():
 p=ROOT/'docs/evidence/p12a_passband_followup_20260926';sample=p/'SAME_GALAXY_SAMPLE.npz';a=np.load(sample)['values'];rows=[]
 for lo,hi in zip([.15,.25,.35,.45],[.25,.35,.45,.55]):
  b=a[(a[:,0]>=lo)&(a[:,0]<hi)];shift=-.8*(b[:,0]-.1);res=b[:,1]+shift
  rows.append(dict(z=[lo,hi],n=len(b),median_sdss_minus_decam_r=float(np.median(b[:,1])),median_existing_term=float(np.median(shift)),median_conditional_residual=float(np.median(res)),residual_p16_p84=np.quantile(res,[.16,.84]).tolist()))
 files=[Path('/global/u2/d/dkololgi/prepare_mocks_Y3_bright.py'),Path('/global/common/software/desi/perlmutter/desiconda/current/code/LSS/main/scripts/mock_tools/prepare_mocks_Y3_bright.py'),ROOT/'docs/evidence/p12a_photometry_recipe_20260926/source/cut_sky_evolution.py',ROOT/'docs/evidence/p12a_photometry_recipe_20260926/source/luminosity_function.py']
 sources=[]
 for i,f in enumerate(files):
  s=f.read_text();name=f'{i:02d}_{f.name}.txt';(OUT/name).write_text(s);sources.append(dict(path=str(f),snapshot=name,sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
 result=dict(rows=rows,assumption='ONLY IF -0.8(z-.1) is an apparent-magnitude adjustment to an otherwise SDSS-like magnitude at fixed physical M, it can offset some/all of the same-SED SDSS-DECam r difference. If R_MAG_ABS is evolution-corrected instead, this arithmetic is not a photometric adjustment. Neither hypothesis is established.',sources=sources,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),sample_sha256=hashlib.sha256(sample.read_bytes()).hexdigest(),no_repairs=True)
 (OUT/'HYPOTHESIS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
