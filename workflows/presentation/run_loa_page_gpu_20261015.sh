#!/bin/bash
# Bounded compute-node payload. The owning salloc releases on exit; no retries.
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export CUBLAS_WORKSPACE_CONFIG=:4096:8
cd /global/u2/d/dkololgi/GraphWeb_DESI
pilot_python=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
"$pilot_python" -u workflows/presentation/loa_model_pilot_20261015.py infer
"$pilot_python" -u workflows/presentation/export_loa_model_page_20261015.py
