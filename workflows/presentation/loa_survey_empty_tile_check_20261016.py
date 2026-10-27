"""One real zero-retained-galaxy tile; run inside the approved GPU preflight.

No job submissions, synthetic truth, production sweep or retries.
"""
import os
import subprocess
import sys
import time
import numpy as np
import loa_survey_runner_20261016 as runner


def main():
    runner.require_compute()
    root=runner.DEFAULT_ROOT; index=runner.read(root/'INDEX.json')
    case=next(r for r in index['cases'] if r['retained_rows']==0)
    started=time.monotonic()
    for action in ('prepare','infer'):
        command=[sys.executable,'-u',runner.__file__,action,'--root',str(root),
                 '--tiles',case['id'],'--stop','1','--seconds','180']
        subprocess.run(command,check=True)
    records=[]
    for seed in runner.SEEDS:
        coarse_path=root/'tiles'/case['id']/f'coarse_{seed}/000'
        coarse=runner.read(coarse_path/'COMPLETE.json')
        runner.verified_bundle(coarse_path,coarse['binding'])
        for architecture,prediction in runner.FAMILIES:
            model=f'{architecture}_{prediction}_{seed}'
            path=root/'tiles'/case['id']/model/'000'
            record=runner.verified_bundle(path,runner.digest([coarse['binding'],model,0,coarse['sha256']]))
            with np.load(path/'arrays.npz') as arrays:
                if arrays['galaxy_eigen'].shape!=(0,3) or arrays['galaxy_rho'].shape!=(0,):
                    raise ValueError('Empty galaxy output has incorrect shape')
                if arrays['rho'].shape!=(32,16,16) or not np.isfinite(arrays['tensor']).all():
                    raise ValueError('Empty-galaxy tile lost its field output')
            records.append(record)
    result=dict(passed=True,tile=case['id'],parent_rows=case['rows'],retained_rows=0,
                models=len(records),records=records,seconds=time.monotonic()-started,
                job=os.environ['SLURM_JOB_ID'],production_launched=False)
    runner.atomic_json(root/'EMPTY_TILE_VALIDATION.json',result)
    print('EMPTY_TILE_VALIDATION_PASSED',flush=True)


if __name__=='__main__':
    main()
