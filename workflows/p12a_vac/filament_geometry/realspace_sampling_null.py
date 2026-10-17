"""Finite-count covariance reference: isotropic uniform spheres, not mock truth."""
import json,os,sys
from pathlib import Path
import numpy as np
from compare_offslab import shape
from prepare_parent import digest
from mock_truth import R
root=Path(sys.argv[1]);assert os.environ.get('SLURM_JOB_ID');report=json.loads((root/'COMPARISON.json').read_text());rng=np.random.default_rng(171020);rows=[]
for r in report['regions']:
 n=r['observed_count_target'];x=rng.normal(size=(512,n,3));x/=np.linalg.norm(x,axis=2)[...,None];x*=R*rng.random((512,n,1))**(1/3)
 s=np.array([shape(a)[0] for a in x]);assert np.isfinite(s).all();np.testing.assert_allclose(s.sum(1),1,atol=1e-12)
 rows.append({'key':r['key'],'n':n,'draws':512,'LPS_median':np.median(s,axis=0).tolist(),'dominant_counts':np.bincount(s.argmax(1),minlength=3).tolist(),'L_95_range':np.quantile(s[:,0],[.025,.975]).tolist()})
(root/'ISOTROPIC_NULL.json').write_text(json.dumps({'seed':171020,'rows':rows,'script_sha256':digest(__file__),'interpretation':'Post-analysis finite-count reference. Covariance line dominance is not a significance test. Uniform real-space sphere without galaxy clustering or DESI selection.'},indent=2)+'\n');print(json.dumps(rows),flush=True)
