#!/bin/bash
# Run as the foreground command of an approved CPU salloc allocation.
set -e
REPO=/global/u2/d/dkololgi/GraphWeb_DESI
OUT=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/v1_forward_ph000_20260928
BASE=$OUT/SecondGenMocks/AbacusSummitBGS_v2
source /global/common/software/desi/desi_environment.sh main >/dev/null 2>&1
module load LSS/main
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
python "$REPO/workflows/p12a_vac/prepare_v1_forward.py" --mockver ab_secondgen_cosmosim --mockpath /global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky --realmin 0 --realmax 1 --prog bright --base_output "$OUT" --apply_mask y --downsampling n --isProduction n --rbandcut 19.5 > "$OUT/prepare_foreground.log" 2>&1
(
 unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
 export PYTHONNOUSERSITE=1
 /pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python "$REPO/workflows/p12a_vac/check_v1_preparation.py" "$BASE/forFA0.fits"
) > "$OUT/validation.log" 2>&1
python "$REPO/workflows/p12a_vac/initialize_v1_forward.py" "$BASE/forFA0.fits" "$BASE/altmtl0" BRIGHT > "$OUT/initialize_foreground.log" 2>&1
python "$REPO/workflows/p12a_vac/run_v1_altmtl.py" --root "$BASE" --smoke > "$OUT/smoke_foreground.log" 2>&1
