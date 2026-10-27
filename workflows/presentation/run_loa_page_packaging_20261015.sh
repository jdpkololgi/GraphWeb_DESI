#!/bin/bash
# Explicitly approved additional 10-minute CPU packaging check; no inference.
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /global/u2/d/dkololgi/GraphWeb_DESI
exec /pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python -u workflows/presentation/export_loa_model_page_20261015.py
