"""Summarise observed read-only closure evidence without modifying samples."""
from pathlib import Path
import json,hashlib
import numpy as np
P=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_closure_20260926'
def main():
 x=json.loads((P/'PRODUCT_CROSSWALK.json').read_text());r=json.loads((P/'MOCK_RETENTION.json').read_text());t=json.loads((P/'TARGET_REPLAY.json').read_text());co=lambda a:np.array(a).reshape(2,4,10).sum(-1)
 counts={k:co(v['counts']) for k,v in x['products'].items()};stage={k:co(v) for k,v in r['stage'].items()}
 assert t['rows']==t['matched']==t['matched_pass']==5436413;assert t['replay_stored_mismatch']==0
 for k in ['kibo','loa']:assert np.array_equal(counts['mock_'+k],stage[k+'_zwarn0'])
 ratios={a+'_over_'+b:(counts[a]/counts[b]).tolist() for a,b in [('real_loa','real_kibo'),('mock_loa','mock_kibo'),('real_loa','mock_loa'),('real_kibo','mock_kibo')]}
 for key in ['kibo','loa']:ratios[key+'_assignment_retention']=(stage[key+'_assigned']/stage[key+'_all']).tolist();ratios[key+'_all_over_forFA']=(stage[key+'_all']/stage['forFA_bright']).tolist()
 out=dict(counts={k:v.tolist() for k,v in counts.items()},stage={k:v.tolist() for k,v in stage.items()},ratios=ratios,checks_pass=True,scope='Diagnostics only; generator receipt and completeness closure still incomplete.',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest());(P/'SUMMARY.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
