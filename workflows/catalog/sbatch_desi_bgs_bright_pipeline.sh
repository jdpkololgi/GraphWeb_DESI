#!/bin/bash
#SBATCH -J desi_bgs_bright
#SBATCH -A desi
#SBATCH -q regular
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 128
#SBATCH -C cpu
#SBATCH -o /pscratch/sd/d/dkololgi/graphweb_desi/logs/desi_bgs_bright_pipeline/%x_%j.out
#SBATCH -e /pscratch/sd/d/dkololgi/graphweb_desi/logs/desi_bgs_bright_pipeline/%x_%j.err
#
# Run one pipeline stage on a compute node (never the login node).
#
#   STEP=2  (default)  Gudhi full graph — CPU, 12h
#   STEP=2b            cuGraph GNN features — submit with GPU overrides (see below)
#   STEP=3             Wedge subset — CPU, 2h typical
#
# Examples:
#   J2=$(sbatch --parsable --export=ALL,STEP=2 /global/homes/d/dkololgi/GraphWeb_DESI/workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh)
#   J2b=$(sbatch --parsable -d afterok:${J2} -A desi -q regular -t 06:00:00 -N 1 -n 1 -c 32 -C gpu --gpus-per-node=1 \
#     --export=ALL,STEP=2b /global/homes/d/dkololgi/GraphWeb_DESI/workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh)
#   sbatch -d afterok:${J2b} --export=ALL,STEP=3 -t 04:00:00 \
#     /global/homes/d/dkololgi/GraphWeb_DESI/workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh
#
# Chain submit from login (lightweight):
#   SUBMIT_CHAIN=1 bash /global/homes/d/dkololgi/GraphWeb_DESI/workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh

set -euo pipefail

if [[ "${SUBMIT_CHAIN:-0}" == "1" ]]; then
  SCRIPT="/global/homes/d/dkololgi/GraphWeb_DESI/workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh"
  LOGDIR="/pscratch/sd/d/dkololgi/graphweb_desi/logs/desi_bgs_bright_pipeline"
  mkdir -p "${LOGDIR}"
  J2=$(sbatch --parsable --export=ALL,STEP=2 "${SCRIPT}")
  echo "Submitted STEP=2 job ${J2}"
  J2B=$(sbatch --parsable -d "afterok:${J2}" -J desi_bgs_cugraph -A desi -q regular -t 06:00:00 \
    -N 1 -n 1 -c 32 -C gpu --gpus-per-node=1 \
    -o "${LOGDIR}/desi_bgs_cugraph_%j.out" -e "${LOGDIR}/desi_bgs_cugraph_%j.err" \
    --export=ALL,STEP=2b "${SCRIPT}") || \
  J2B=$(sbatch --parsable -d "afterok:${J2}" -J desi_bgs_cugraph -A ntrain6 -q regular -t 06:00:00 \
    -N 1 -n 1 -c 32 -C gpu --gpus-per-node=1 \
    -o "${LOGDIR}/desi_bgs_cugraph_%j.out" -e "${LOGDIR}/desi_bgs_cugraph_%j.err" \
    --export=ALL,STEP=2b "${SCRIPT}")
  echo "Submitted STEP=2b job ${J2B} (afterok ${J2})"
  J3=$(sbatch --parsable -d "afterok:${J2B}" -J desi_bgs_wedge -A desi -q regular -t 04:00:00 \
    -N 1 -n 1 -c 64 -C cpu \
    -o "${LOGDIR}/desi_bgs_wedge_%j.out" -e "${LOGDIR}/desi_bgs_wedge_%j.err" \
    --export=ALL,STEP=3 "${SCRIPT}") || \
  J3=$(sbatch --parsable -d "afterok:${J2B}" -J desi_bgs_wedge -A ntrain6 -q regular -t 04:00:00 \
    -N 1 -n 1 -c 64 -C cpu \
    -o "${LOGDIR}/desi_bgs_wedge_%j.out" -e "${LOGDIR}/desi_bgs_wedge_%j.err" \
    --export=ALL,STEP=3 "${SCRIPT}")
  echo "Submitted STEP=3 job ${J3} (afterok ${J2B})"
  echo "LOGDIR=${LOGDIR}"
  exit 0
fi

STEP="${STEP:-2}"

REPO="/global/homes/d/dkololgi/GraphWeb_DESI"
SCRATCH="/pscratch/sd/d/dkololgi/graphweb_desi"
CATALOG="${SCRATCH}/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits"
GRAPH_DIR="${SCRATCH}/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright"
WEDGE_DIR="${SCRATCH}/outputs/desi_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_bright_from_fullgraph"
WEDGE_PREFIX="desi_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_bright"

PY_COSMIC="/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python"
PY_RAPIDS="/pscratch/sd/d/dkololgi/conda/envs/rapids-gnn/bin/python"

module purge 2>/dev/null || true
unset PYTHONPATH PYTHONHOME
export PYTHONNOUSERSITE=1

mkdir -p "${GRAPH_DIR}" "${WEDGE_DIR}"

echo "=== desi_bgs_bright pipeline STEP=${STEP} ==="
echo "host=$(hostname) job=${SLURM_JOB_ID:-none} step=${STEP}"
date

run_step2() {
  srun -n 1 -c "${SLURM_CPUS_PER_TASK:-128}" --cpu-bind=cores \
    "${PY_COSMIC}" "${REPO}/workflows/graph_construction/build_desi_bgs_gudhi_graph.py" \
    --catalog-path "${CATALOG}" \
    --out-dir "${GRAPH_DIR}" \
    --alpha-sq inf \
    --seed 42 \
    --chunk-rows 2000000 \
    --output-prefix desi_delaunay \
    --heartbeat-seconds 60
}

run_step2b() {
  CU_META="${GRAPH_DIR}/graph_metadata.json"
  [[ -f "${CU_META}" ]] || CU_META="${GRAPH_DIR}/desi_delaunay_metadata.json"
  srun -n 1 -c "${SLURM_CPUS_PER_TASK:-32}" --gpus-per-task=1 \
    "${PY_RAPIDS}" "${REPO}/workflows/graph_construction/desi_graph_features_cugraph.py" \
    --metadata-path "${CU_META}" \
    --out-dir "${GRAPH_DIR}" \
    --out-prefix desi_bgs_cugraph
}

run_step3() {
  WEDGE_META="${GRAPH_DIR}/desi_delaunay_metadata.json"
  [[ -f "${WEDGE_META}" ]] || WEDGE_META="${GRAPH_DIR}/graph_metadata.json"
  srun -n 1 -c "${SLURM_CPUS_PER_TASK:-64}" --cpu-bind=cores \
    "${PY_COSMIC}" "${REPO}/workflows/graph_construction/subset_desi_graph_wedge.py" \
    --graph-metadata "${WEDGE_META}" \
    --catalog-path "${CATALOG}" \
    --out-dir "${WEDGE_DIR}" \
    --out-prefix "${WEDGE_PREFIX}" \
    --ra-min 120 --ra-max 140 \
    --dec-min 16.5 --dec-max 26.7 \
    --z-min 0.25 --z-max 0.3 \
    --parent-gnn-arrays "${GRAPH_DIR}/desi_bgs_cugraph_gnn_arrays.npz" \
    --parent-gnn-metadata "${GRAPH_DIR}/desi_bgs_cugraph_gnn_metadata.json"
}

case "${STEP}" in
  2)   run_step2 ;;
  2b)  run_step2b ;;
  3)   run_step3 ;;
  *)
    echo "Unknown STEP=${STEP} (use 2, 2b, or 3)" >&2
    exit 2
    ;;
esac

date
echo "=== STEP=${STEP} finished ==="
