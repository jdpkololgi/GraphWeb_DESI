#!/bin/bash
# Submit only after the preview array completes, or with its exact afterok ID.
# This audit does not approve or submit the 32-draw extension.
#SBATCH --job-name=loa6-preview-audit
#SBATCH --account=desi
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --licenses=scratch
#SBATCH --no-requeue
#SBATCH --output=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1/logs/preview_audit_%j.log
#SBATCH --error=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1/logs/preview_audit_%j.err
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
cd /global/u2/d/dkololgi/GraphWeb_DESI/workflows/presentation
srun --nodes=1 --ntasks=1 --cpus-per-task=64 --cpu-bind=cores \
  /pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python -u audit_loa_survey_preview_20261016.py
