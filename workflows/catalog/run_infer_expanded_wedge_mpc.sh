#!/bin/bash
# Jraph inference on Mpc-parity DESI expanded wedge (after rerun_desi_bgs_bright_graph_wedge_mpc.sh).
set -euo pipefail

REPO="/global/homes/d/dkololgi/GraphWeb_DESI"
# shellcheck source=/dev/null
source "${REPO}/workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt"

export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
export ILLUSTRIS_ROOT="${ILLUSTRIS_ROOT:-/global/homes/d/dkololgi/TNG/Illustris}"
export JAX_PLATFORMS=cpu

"${PY_COSMIC}" "${REPO}/workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py" \
  --abacus-run-dir "${ABACUS_RUN_DIR}" \
  --calibration-cache "${CALIBRATION_CACHE}" \
  --abacus-gnn-arrays "${ABACUS_WEDGE_GNN_ARRAYS}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/inference_outputs \
  --run-name "${INFER_RUN_NAME}" \
  --latent-size "${LATENT_SIZE}" \
  --num-heads "${NUM_HEADS}" \
  --num-passes "${NUM_PASSES}" \
  --dropout "${DROPOUT}" \
  --seed "${SEED}" \
  --lambda-threshold "${LAMBDA_THRESHOLD}"

echo "Wrote: ${INFER_DIR}"
