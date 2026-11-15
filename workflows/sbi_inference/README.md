# FlowJAX/SBI DESI wedge inference

This directory contains the uncertainty-aware DESI wedge inference path. It is
the SBI/NPE companion to the Jraph regression workflow: both consume the
expanded Mpc-parity DESI wedge graph arrays, but FlowJAX returns posterior
samples, posterior widths, and T-Web class probabilities for each galaxy.

## When to use this path

Use `infer_desi_wedge_flowjax.py` when the analysis needs:

- per-galaxy posterior means for `(lambda_1, lambda_2, lambda_3)`,
- posterior widths as uncertainty diagnostics,
- class probabilities for void/wall/filament/cluster,
- truth-free DESI checks against Abacus training/reference distributions.

Use `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py` instead
when you only need the deterministic Jraph regression point estimate.

## Required inputs

The DESI-side inputs are the expanded wedge products listed in
`workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt`:

- `${DESI_GNN_ARRAYS}`
- `${DESI_GNN_METADATA}`
- `${DESI_GLOBAL_NODE_IDS}`
- `${DESI_WEDGE_CATALOG_NPZ}`

The inference code also needs Abacus/Illustris training artifacts:

- `ILLUSTRIS_ROOT`: sibling Illustris repo containing the FlowJAX/GNN code.
- `--model-path`: FlowJAX SBI checkpoint.
- `--calibration-cache`: matching training cache with `node_feature_scaler` and
  Abacus eigenvalues.
- `--abacus-gnn-arrays`: path1 fiberassign Abacus wedge `.npz` used to fit the
  edge scaler for the FlowJAX training graph.

Do not mix model/cache variants. If the model was trained with
scale-invariant features, pass `--scale-invariant-features` and use the matching
scale-invariant cache.

## Inference command

Run on a GPU node with the same environment used for JAX/FlowJAX work:

```bash
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

export ILLUSTRIS_ROOT=/global/homes/d/dkololgi/TNG/Illustris
export FLOWJAX_MODEL_PATH=/path/to/flowjax_sbi_model_seed_42_best.pkl
export FLOWJAX_CACHE=/path/to/processed_jraph_data_scaled_linear_eig.pkl
export PATH1_EDGE_SCALER_NPZ=/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz

python workflows/sbi_inference/infer_desi_wedge_flowjax.py \
  --model-path "${FLOWJAX_MODEL_PATH}" \
  --calibration-cache "${FLOWJAX_CACHE}" \
  --abacus-gnn-arrays "${PATH1_EDGE_SCALER_NPZ}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --num-posterior-samples 128 \
  --lambda-threshold 0.2 \
  --scale-invariant-features \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs \
  --run-name desi_wedge_flowjax_linear_si
```

The script writes:

- `desi_wedge_flowjax_preds.npz`
  - `global_node_id`, `ra`, `dec`, `z`
  - `lambda_mean[N,3]`, `lambda_std[N,3]`
  - `classprob[N,4]`, `hard_class[N]`, `p_exceed[N,3]`
  - GNN embeddings and a subset of posterior samples for plotting
- `summary.json`
  - model/cache provenance,
  - node/edge scaler diagnostics,
  - T-Web class fractions,
  - ordering and class-probability consistency checks.

## Diagnostics and figures

Generate DESI figures after inference:

```bash
PRED_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear_si

python workflows/sbi_inference/plot_desi_wedge_flowjax.py \
  --preds-npz "${PRED_DIR}/desi_wedge_flowjax_preds.npz" \
  --summary-json "${PRED_DIR}/summary.json" \
  --calibration-cache "${FLOWJAX_CACHE}" \
  --abacus-self-npz /path/to/abacus_self_flowjax_preds.npz \
  --output-dir "${PRED_DIR}/figures"
```

Useful plot flags:

- `--abacus-classprob-npz`: adds Abacus NPE class fractions.
- `--abacus-self-npz`: enables Abacus-vs-DESI eigenvalue and embedding overlays.
- `--calibration-cache`: supplies Abacus truth eigenvalues when a self-inference
  `.npz` is not available.
- `--no-regression`: omits deterministic regression bars from the class-fraction
  panel.
- `--only-class-fractions`: regenerates only `class_fractions_comparison.png`.

Build an Abacus self-reference file for the overlays:

```bash
python workflows/sbi_inference/infer_abacus_self_flowjax.py \
  --model-path "${FLOWJAX_MODEL_PATH}" \
  --calibration-cache "${FLOWJAX_CACHE}" \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear_si
```

## Property/environment closure

`build_desi_wedge_property_join.py` joins the NPE predictions to LOA FastSpecFit
properties by recovering `TARGETID` from `global_node_id`. The
`global_node_id` values index the exact source BGS FITS catalog used to build the
graph, so `--source-catalog` must match the wedge construction input.

Duplicate TARGETIDs from the Galactic N/S graph split are averaged across their
posterior products, then the hard class is recomputed from the averaged class
probabilities. The writer prefers Parquet and falls back to CSV if the Parquet
engine is unavailable.

```bash
python workflows/sbi_inference/build_desi_wedge_property_join.py \
  --preds "${PRED_DIR}/desi_wedge_flowjax_preds.npz"

python workflows/sbi_inference/plot_property_environment_closure.py \
  --table "${PRED_DIR}/desi_wedge_env_props.parquet"
```

## Operational constraints and pitfalls

- DESI graph metadata must use `coordinate_units: "mpc"`. Mpc/h or unknown
  metadata aborts because Abacus training parity would be broken.
- The DESI node-feature columns must be exactly `Degree`, `Clustering`,
  `Density`, `Neigh Density`, `I_eig1`, `I_eig2`, `I_eig3`.
- The path1 fiberassign wedge is required for `--abacus-gnn-arrays`; same-sized
  regression wedges can silently pass shape checks, so the script asserts scaler
  constants in non-scale-invariant mode.
- `--edge-domain-adapt` and `--node-domain-adapt` are diagnostic domain-shift
  experiments. Do not treat them as production defaults without recording the
  science decision in `SCIENCE_LOG.md`.
- Posterior samples are sorted into ascending eigenvalue order by default. Use
  `--no-sort` only to inspect raw flow ordering violations.
- `plot_desi_wedge_flowjax.py` writes
  `eigenvalue_distributions_3way.png`; older comments and print output may refer
  to this diagnostic as `domain_shift_overlay.png`.
- `scripts/sync_figures_to_canonical.sh` mirrors top-level figure files from a
  run directory into `$GRAPHWEB_CANONICAL_FIGURE_DIR/<run-name>/`.
