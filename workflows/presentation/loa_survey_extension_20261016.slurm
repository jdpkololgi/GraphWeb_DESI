#!/bin/bash
# Prepared only. Do not submit until the complete four-draw preview is reviewed.
#SBATCH --job-name=loa6-extend32
#SBATCH --account=desi_g
#SBATCH --qos=shared
#SBATCH --constraint=gpu&hbm80g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus=1
#SBATCH --time=10:00:00
#SBATCH --array=0-102%8
#SBATCH --licenses=scratch
#SBATCH --no-requeue
#SBATCH --output=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1/logs/extension_%A_%a.log
#SBATCH --error=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1/logs/extension_%A_%a.err
set -euo pipefail
srun --nodes=1 --ntasks=1 --cpus-per-task=32 --gpus=1 --cpu-bind=cores \
  bash /global/u2/d/dkololgi/GraphWeb_DESI/workflows/presentation/run_loa_survey_production_20261016.sh \
  shard --shard "$SLURM_ARRAY_TASK_ID" --draws 32
