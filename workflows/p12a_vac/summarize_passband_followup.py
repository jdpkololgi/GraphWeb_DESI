"""Summarize sensitivity subsets and validate saved diagnostic provenance."""
from pathlib import Path
import json,hashlib,ast
import numpy as np
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_passband_followup_20260926'
def main():
 x=np.load(OUT/'SAME_GALAXY_SAMPLE.npz');a=x['values'];out=[]
 for lo,hi in zip([.15,.25,.35,.45],[.25,.35,.45,.55]):
  b=a[(a[:,0]>=lo)&(a[:,0]<hi)];row={'z':[lo,hi]}
  for key,k in [('all',np.ones(len(b),bool)),('colour_residual_lt003',abs(b[:,5])<.03),('rchi2phot_lt10',b[:,-1]<10)]:
   c=b[k];row[key]=dict(n=len(c),median_sdss_minus_decam_r=float(np.median(c[:,1])),median_sdss_minus_decam_gr=float(np.median(c[:,2]))) if len(c) else dict(n=0)
  out.append(row)
 (OUT/'SED_ROBUSTNESS.json').write_text(json.dumps(out,indent=2)+'\n')
 checks=[]
 for fn,sn in [('SAME_GALAXY','same_galaxy_passbands'),('SDSS_IMAGING_PAIRS','sdss_imaging_pairs'),('CATALOGUE_PARITY','passband_catalogue_parity')]:
  d=json.loads((OUT/(fn+'.json')).read_text());p=Path(__file__).with_name(sn+'.py');ast.parse(p.read_text());assert d['script_sha256']==hashlib.sha256(p.read_bytes()).hexdigest();checks.append(fn+' source hash and syntax')
 assert len(np.unique(x['targetid']))==len(x['targetid']);assert np.isfinite(a).all();assert np.max(abs(a[:,3:5]))<.00002;checks+=['unique target IDs, finite synthesis, model-colour replay below 0.00002 mag']
 d=json.loads((OUT/'SDSS_IMAGING_PAIRS.json').read_text());assert len(set(r['targetid'] for r in d['records']))==len(d['records']);assert all(r['sep_arcsec']<1 for r in d['records']);checks+=['unique imaging matches below 1 arcsec']
 parity=json.loads((OUT/'CATALOGUE_PARITY.json').read_text());assert max(parity['correct_flux_relative_residual'].values())<1e-5;assert parity['correct_fsf_minus_loa_dered_gr']['max_abs']<.001;checks+=['corrected FastSpecFit flux matches Loa divided once by FastSpecFit MW transmission; colour difference below .001mag']
 (OUT/'VALIDATION.json').write_text(json.dumps(dict(passed=True,checks=checks,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='VALIDATION.json'}),indent=2)+'\n')
 print('\n'.join(checks))
if __name__=='__main__':main()
