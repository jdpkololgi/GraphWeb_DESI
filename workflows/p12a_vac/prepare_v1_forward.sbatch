#!/bin/bash
#SBATCH --job-name=v1-ph000-prepare
#SBATCH --account=desi
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --licenses=scratch
#SBATCH --output=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/v1_forward_ph000_20260928/batch-%j.log
set -eo pipefail
REPO=/global/u2/d/dkololgi/GraphWeb_DESI
OUT=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/v1_forward_ph000_20260928
cd "$REPO"
sha256sum -c docs/evidence/p12a_joint_colour_luminosity_20260928/FORWARD_SOURCE.sha256
source /global/common/software/desi/desi_environment.sh main >/dev/null 2>&1
module load LSS/main
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
python - <<'PY'
import LSS,subprocess
from pathlib import Path
p=Path(LSS.__file__).resolve().parents[2]
a=subprocess.check_output(['git','-c','safe.directory='+str(p),'-C',str(p),'rev-parse','HEAD'],text=True).strip()
assert a=='d942b990860e016e7558c740966fea3e5ce7e7f3',a
PY
srun --ntasks=1 --cpus-per-task=64 python workflows/p12a_vac/prepare_v1_forward.py --mockver ab_secondgen_cosmosim --mockpath /global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky --realmin 0 --realmax 1 --prog bright --base_output "$OUT" --apply_mask y --downsampling n --isProduction n --rbandcut 19.5
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
srun --ntasks=1 --cpus-per-task=2 /pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python workflows/p12a_vac/check_v1_preparation.py "$OUT/SecondGenMocks/AbacusSummitBGS_v2/forFA0.fits"
