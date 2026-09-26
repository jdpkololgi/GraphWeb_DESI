#!/bin/bash
#SBATCH --account=desi_g
#SBATCH --constraint=gpu&hbm80g
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus=1
#SBATCH --mem=48G
#SBATCH --time=00:45:00
#SBATCH --licenses=scratch
set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
export ILLUSTRIS_ROOT=/pscratch/sd/d/dkololgi/abacus/p10_multiphase/p12a_halo48_candidate_20260924_v1/source
VAC_RUN_ROOT=${1:?production root required}
VAC_STAGE=${2:-worker}
VAC_PYTHON=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
if [[ "$VAC_STAGE" == worker ]]; then
    srun --mpi=none --resv-ports=0 --ntasks=1 --cpus-per-task=32 --gpus=1 "$VAC_PYTHON" "$VAC_RUN_ROOT/source/workflows/p12a_vac/production.py" worker --root "$VAC_RUN_ROOT" --part "${SLURM_ARRAY_TASK_ID:?array index required}"
else
    srun --mpi=none --resv-ports=0 --ntasks=1 --cpus-per-task=32 --gpus=1 "$VAC_PYTHON" "$VAC_RUN_ROOT/source/workflows/p12a_vac/production.py" merge --root "$VAC_RUN_ROOT"
fi
