"""Measured historical observation response applied as a transfer diagnostic only."""
from pathlib import Path
import json,hashlib
import numpy as np
R=Path(__file__).resolve().parents[2]
def main():
 p=R/'docs/evidence/p12a_alignment_execution_20260926/parent'
 v1=np.load(p/'v1.npz')['z'];old=np.load(p/'v0.1.npz')['z']
 prep=json.load(open('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_mock_test_ph006_20260926_v1/PREPARED.json'))
 stage=json.load(open(R/'docs/evidence/p12a_closure_20260926/MOCK_RETENTION.json'))['stage']
 raw=np.array(prep['raw_counts']);response=np.array(stage['loa_zwarn0'])/raw
 assert raw.shape==v1.shape==(2,40)
 assert np.all((response>=0)&(response<=1))
 shell=lambda x:np.asarray(x).reshape(2,4,10).sum(-1)
 targets={k:shell(v) for k,v in prep['stages'].items()}
 transferred=shell(v1*response)
 out={'qualification':False,'method':'ph006 old-parent observed/raw ratio per cap and dz=.01 applied to ph000 v1. Diagnostic assumption, NOT v1 fibreassign or a validated quality model.', 'target_counts':{k:v.tolist() for k,v in targets.items()},'raw_v1_counts':shell(v1).tolist(),'transferred_v1_counts':transferred.tolist(),'old_response_shell':(shell(np.array(stage['loa_zwarn0']))/shell(raw)).tolist(),'ratios':{k:(transferred/v).tolist() for k,v in targets.items()},'missing_v1_observables':['DELTACHI2','SPECTYPE','ZWARN','fibre magnitude','imaging/targeting nuisance observables','actual Loa assignment and HPmap/veto output'], 'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
 (R/'docs/evidence/p12a_photometry_swap_20260928/OBSERVATION_SCREEN.json').write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
