#!/bin/bash
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export CUBLAS_WORKSPACE_CONFIG=:4096:8
cd /global/u2/d/dkololgi/GraphWeb_DESI
exec /pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python -u workflows/presentation/loa_survey_production_20261016.py "$@"
