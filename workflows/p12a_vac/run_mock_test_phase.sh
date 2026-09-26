#!/bin/bash
# Run inside an approved CPU allocation. Each stage refuses existing outputs.
set -euo pipefail
cd /global/u2/d/dkololgi/GraphWeb_DESI
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1 PYTHONUNBUFFERED=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export DESIMODEL=/global/common/software/desi/perlmutter/desiconda/current/code/desimodel/main
TEST_PY=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
TEST_ROOT=${2:?provide isolated scratch root}
case ${1:?provide baseline, prepare, sweep or verify} in
 baseline)
  "$TEST_PY" -c 'import numpy as np,runpy; np.random.seed(62026); runpy.run_path("docs/evidence/p12a_nz_magnitude_20260926/upstream_prepare_mocks_Y3_bright.py.txt",run_name="__main__")' \
   --mockver ab_secondgen_cosmosim --realmin 6 --realmax 7 --prog bright \
   --isProduction n --base_output "$TEST_ROOT/upstream_baseline" \
   --apply_mask n --downsampling n --rbandcut 19.5 ;;
 prepare|sweep) "$TEST_PY" workflows/p12a_vac/mock_selection_test.py "$1" --root "$TEST_ROOT" ;;
 verify) "$TEST_PY" workflows/p12a_vac/verify_mock_test_baseline.py --root "$TEST_ROOT" ;;
 *) exit 2 ;;
esac
