#!/bin/bash
#SBATCH --job-name=loa_field_properties
#SBATCH --account=desi
#SBATCH --qos=regular
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH --licenses=scratch
#SBATCH --output=/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1_properties_%j.log

set -euo pipefail
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=/tmp/loa-properties-mpl
ATLAS_ROOT=/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1
ANALYSIS_SOURCE=${1:?Pass the immutable GraphWeb source snapshot directory}
PYTHON=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
cd "$ANALYSIS_SOURCE"
"$PYTHON" - "$ATLAS_ROOT" <<'PY'
import hashlib,json,sys
from pathlib import Path
root=Path(sys.argv[1]);record=json.loads((root/'ATLAS_COMPLETE.json').read_text())
assert record['passed'] and record['rows']==5436413 and record['tiles']==2472 and record['draws']==8
for name,sha in json.loads(Path('SOURCE_SHA256.json').read_text()).items():
    assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha, name
print('Complete atlas and frozen analysis source verified',flush=True)
PY
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$PYTHON" -u workflows/p12a_vac/field_property_environment.py \
  --root "$ATLAS_ROOT" --out "$ATLAS_ROOT/property_comparison_v1"
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$PYTHON" -u workflows/p12a_vac/loa_field_viewer.py \
  --root "$ATLAS_ROOT" --out "$ATLAS_ROOT/viewer_full_v1" --surface-library "$ATLAS_ROOT/viewer_dependencies"
export PLAYWRIGHT_BROWSERS_PATH="$ATLAS_ROOT/browser_runtime"
srun --ntasks=1 --cpus-per-task=32 --cpu-bind=cores "$PYTHON" -u workflows/p12a_vac/check_loa_viewer.py \
  --html "$ATLAS_ROOT/viewer_full_v1/Loa_3D.html" --out "$ATLAS_ROOT/viewer_full_v1/qa" \
  --playwright-path "$ATLAS_ROOT/playwright_runtime"
echo 'Property comparison and full local viewer built; inspect numerical results and screenshots before scientific handoff.'
