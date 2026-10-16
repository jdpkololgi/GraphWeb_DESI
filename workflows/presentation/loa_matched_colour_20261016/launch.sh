#!/usr/bin/env bash
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python /global/homes/d/dkololgi/GraphWeb_DESI/workflows/presentation/loa_matched_colour_20261016/run.py
