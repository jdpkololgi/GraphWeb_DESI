#!/usr/bin/env bash
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
REST_SCRIPT_DIR=/global/homes/d/dkololgi/GraphWeb_DESI/workflows/presentation/loa_rest_colour_20261016
REST_PYTHON=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
"$REST_PYTHON" "$REST_SCRIPT_DIR/run.py" --stage all
"$REST_PYTHON" "$REST_SCRIPT_DIR/render.py"
