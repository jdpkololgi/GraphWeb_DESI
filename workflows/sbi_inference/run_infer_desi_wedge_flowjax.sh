#!/bin/bash
# FlowJAX NPE (linear) inference on the real DESI LOA wedge. 1x A100, ~5 min.
set -euo pipefail
WD=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph
PRE=$WD/desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc
RUN=/pscratch/sd/d/dkololgi/abacus/sbi_runs/path1_wedge_flowjax_3d_Bcorrected_linear
MODEL=$(ls -t $RUN/flowjax_sbi_model_seed_42_*.pkl | head -1)
CACHE=/pscratch/sd/d/dkololgi/abacus/sbi_caches/path1_flowjax_3d_lineareig/processed_jraph_data_mc1e+09_v2_scaled_3_linear_eig.pkl
PATH1=/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz
LOG=/pscratch/sd/d/dkololgi/logs/desi_flowjax_infer_$(date +%Y%m%d_%H%M%S).log
echo "$LOG" > /tmp/desi_flowjax_infer_log.txt
echo "=== DESI FlowJAX NPE inference salloc 1xA100 1h at $(date) ===" | tee "$LOG"
salloc --nodes=1 --gpus-per-node=1 --cpus-per-task=32 --constraint=gpu \
       --qos=interactive --time=01:00:00 --account=desi_g \
  srun -n 1 bash -lc '
    set -euo pipefail
    unset PYTHONPATH PYTHONHOME; export PYTHONNOUSERSITE=1
    source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate cosmic_env
    export XLA_PYTHON_CLIENT_PREALLOCATE=false XLA_PYTHON_CLIENT_ALLOCATOR=platform
    export ILLUSTRIS_ROOT=/global/homes/d/dkololgi/TNG/Illustris
    cd /global/u2/d/dkololgi/GraphWeb_DESI
    python -c "import jax; print(\"jax devices:\", jax.devices())"
    python -u workflows/sbi_inference/infer_desi_wedge_flowjax.py \
      --model-path '"$MODEL"' \
      --calibration-cache '"$CACHE"' \
      --abacus-gnn-arrays '"$PATH1"' \
      --desi-gnn-arrays '"$PRE"'_gnn_arrays.npz \
      --desi-gnn-metadata '"$PRE"'_gnn_metadata.json \
      --desi-global-node-ids '"$PRE"'_global_node_ids.npy \
      --desi-wedge-catalog-npz '"$PRE"'_wedge_catalog_minimal.npz \
      --num-posterior-samples 128 --lambda-threshold 0.2 \
      --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs \
      --run-name desi_wedge_flowjax_linear
  ' 2>&1 | tee -a "$LOG"
echo "=== DESI FLOWJAX INFER EXITED at $(date) (rc=$?) ===" | tee -a "$LOG"
