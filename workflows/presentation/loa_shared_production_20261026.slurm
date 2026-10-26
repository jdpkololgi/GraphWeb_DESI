#!/bin/bash
# Times/job names/log paths are explicit per-cap submission arguments.
#SBATCH --account=desi_g
#SBATCH --qos=shared
#SBATCH --constraint=gpu&hbm80g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --gpus=1
#SBATCH --licenses=scratch
#SBATCH --no-requeue

set -euo pipefail
snapshot=${1:?frozen source root required}
manifest_sha=${2:?manifest hash required}
cap=${3:?cap required}
seed=${4:?seed required}
case "$cap:$seed" in NGC:17|NGC:29|SGC:17|SGC:29) ;; *) exit 2 ;; esac
case "$snapshot" in /pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/source_production_v1) ;; *) exit 2 ;; esac
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
export ILLUSTRIS_ROOT="$snapshot/Illustris"
export PATH="/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin:$PATH"
cd "$snapshot/GraphWeb_DESI/workflows/presentation"
python_bin=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
srun --nodes=1 --ntasks=1 --cpus-per-task=32 --gpus=1 --cpu-bind=cores \
  "$python_bin" -m unittest -v test_loa_shared_sixfit_20261026 test_loa_shared_field_io_20261026
# Local AST maintenance is ancillary; a filesystem lock failure cannot certify
# or invalidate the numerical tests. Never bypass the lock or invoke an API.
if [[ "$cap:$seed" == "NGC:17" ]]; then
  srun --nodes=1 --ntasks=1 --cpus-per-task=32 --gpus=1 --cpu-bind=cores \
    bash -c 'cd /global/u2/d/dkololgi/GraphWeb_DESI && graphify update . && graphify global add graphify-out/graph.json --as GraphWeb_DESI' \
    || printf '%s\n' 'WARNING: repository graph refresh failed; see job log.'
fi
srun --nodes=1 --ntasks=1 --cpus-per-task=32 --gpus=1 --cpu-bind=cores \
  "$python_bin" -u loa_shared_production_20261026.py --snapshot "$snapshot" \
  --manifest-sha256 "$manifest_sha" --cap "$cap" --seed "$seed"
