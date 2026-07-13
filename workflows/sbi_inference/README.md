# FlowJAX/SBI DESI wedge inference

This directory contains the DESI-side posterior inference and diagnostics for the
Abacus-trained FlowJAX neural posterior estimator (NPE). It is the current
DESI wedge path for per-galaxy ordered tidal-tensor eigenvalue posteriors:

```text
DESI Mpc-parity wedge graph
  -> Abacus-parity node/edge transforms
  -> Illustris GNN encoder + FlowJAX posterior
  -> per-galaxy lambda samples, class probabilities, and closure figures
```

This is separate from:

- `workflows/graph_inference/graph_catalog.py`: all-sky GAT VAC-style class
  prediction.
- `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`: Jraph
  point-regression for eigenvalues on a DESI wedge.

## Scientific intent

The model infers the ordered T-web eigenvalues
`lambda_1 <= lambda_2 <= lambda_3`. Void/wall/filament/cluster classes are
derived from posterior samples using `lambda_threshold` (default `0.2`), not
treated as the fundamental target. DESI has no per-galaxy truth labels, so the
core DESI-side checks are truth-free:

- class-probability count/marginal consistency;
- posterior width and entropy maps;
- Abacus-vs-DESI eigenvalue and embedding comparisons;
- galaxy-property closure tests after joining FastSpecFit properties.

Calibration tests that require truth (TARP/SBC and Abacus test-set metrics)
remain in the sibling Illustris training workflow.

## Required inputs and parity constraints

`infer_desi_wedge_flowjax.py` requires a DESI wedge built from the Gudhi/cuGraph
Mpc-parity graph products plus the matching Abacus training artifacts.

Required inputs:

- `--model-path`: FlowJAX model pickle from the Illustris SBI run.
- `--calibration-cache`: training cache with `node_feature_scaler`,
  `eigenvalues_raw`, and target scaler metadata.
- `--abacus-gnn-arrays`: raw path1 fiberassign Abacus wedge GNN arrays used to
  fit the DESI edge scaler.
- `--desi-gnn-arrays`: DESI wedge GNN arrays.
- `--desi-gnn-metadata`: DESI wedge metadata JSON.
- `--desi-global-node-ids`: row indices into the source BGS graph catalog.
- `--desi-wedge-catalog-npz`: minimal RA/Dec/z wedge catalog.

Important invariants enforced by the script:

- DESI coordinate units must be comoving `mpc`; `mpc/h` or unknown metadata
  aborts inference.
- Node features must be exactly:
  `Degree, Clustering, Density, Neigh Density, I_eig1, I_eig2, I_eig3`.
- Edge features are duplicated bidirectionally, then edge length and density
  contrast are log-transformed and standardized with an Abacus path1 wedge
  scaler.
- `ILLUSTRIS_ROOT` must point to the sibling Illustris repo so imports such as
  `shared.graph_net_models` and `workflows.sbi.plot_flowjax_posteriors` resolve
  to the training code.

Use the expanded Mpc wedge manifest as the source of canonical DESI paths:

```bash
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt
```

## Inference

The checked-in launcher captures the validated baseline Perlmutter pattern:

```bash
bash workflows/sbi_inference/run_infer_desi_wedge_flowjax.sh
```

The underlying command shape is:

```bash
export ILLUSTRIS_ROOT=/global/homes/d/dkololgi/TNG/Illustris
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

FLOWJAX_RUN=/pscratch/sd/d/dkololgi/abacus/sbi_runs/path1_wedge_flowjax_3d_Bcorrected_linear
MODEL="${FLOWJAX_RUN}/flowjax_sbi_model_seed_42_YYYYMMDD_HHMMSS.pkl"
CACHE=/pscratch/sd/d/dkololgi/abacus/sbi_caches/path1_flowjax_3d_lineareig/processed_jraph_data_mc1e+09_v2_scaled_3_linear_eig.pkl
PATH1_EDGE_ARRAYS=/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz

python -u workflows/sbi_inference/infer_desi_wedge_flowjax.py \
  --model-path "${MODEL}" \
  --calibration-cache "${CACHE}" \
  --abacus-gnn-arrays "${PATH1_EDGE_ARRAYS}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --num-posterior-samples 128 \
  --lambda-threshold 0.2 \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs \
  --run-name desi_wedge_flowjax_linear
```

For the scale-invariant production run, replace `MODEL` and `CACHE` with the
matching scale-invariant training artifacts, set
`--run-name desi_wedge_flowjax_linear_si`, and add
`--scale-invariant-features`.

Outputs are written to `--output-dir/--run-name/`:

- `desi_wedge_flowjax_preds.npz`: `global_node_id`, sky coordinates,
  `lambda_mean`, `lambda_std`, `classprob`, `hard_class`, `p_exceed`,
  embeddings, and a posterior-sample subset.
- `summary.json`: provenance, scaler statistics, class fractions, ordering
  checks, and Abacus reference fractions.

### Domain-adaptation and production flags

The inference script exposes three transfer-diagnostic modes:

- `--edge-domain-adapt`: re-standardizes DESI scaled edge length and density
  contrast to training-like `N(0,1)`.
- `--node-domain-adapt`: re-standardizes DESI box-cox node features to
  training-like `N(0,1)`.
- `--scale-invariant-features`: uses per-graph-median normalized node and edge
  scale features. Use this only with a matching scale-invariant training cache
  and model.

Treat node/edge domain adaptation as diagnostics unless the science log records
the corresponding run as the selected data product. The scale-invariant mode is
the durable Route A transfer path, but only when model, cache, and DESI transform
all match.

## Abacus self-inference reference

For DESI comparison plots, run the trained NPE on its own Abacus test wedge:

```bash
python workflows/sbi_inference/infer_abacus_self_flowjax.py \
  --model-path "${MODEL}" \
  --calibration-cache "${CACHE}" \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear_si
```

This writes `abacus_self_flowjax_preds.npz`, which supplies Abacus NPE
posterior means, embeddings, and class probabilities for three-way
Abacus-truth / Abacus-NPE / DESI-NPE figures.

## Figures and truth-free diagnostics

After DESI inference, generate the main diagnostic figures on a CPU node:

```bash
RUN_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear_si
ABACUS_SELF=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear_si/abacus_self_flowjax_preds.npz

python workflows/sbi_inference/plot_desi_wedge_flowjax.py \
  --preds-npz "${RUN_DIR}/desi_wedge_flowjax_preds.npz" \
  --summary-json "${RUN_DIR}/summary.json" \
  --calibration-cache "${CACHE}" \
  --abacus-self-npz "${ABACUS_SELF}" \
  --output-dir "${RUN_DIR}" \
  --no-regression
```

Common products include:

- `class_fractions_comparison.png`;
- `class_sky_map.png`;
- `posterior_width_sky_map.png`;
- `eigenvalue_distributions_3way.png`;
- `embedding_pca.png` and optional `embedding_pca_3d.html`;
- optional UMAP and 3D HTML figures when dependencies are available.

For the cluster-deficit transfer story, see:

- `CLUSTER_DEFICIT_FIX_PLAN.md`;
- `plot_edge_scale.py`;
- `plot_cluster_recovery_bars.py`;
- `plot_eigenvalue_corner.py`;
- `plot_eig_dist_buildup.py`.

## Galaxy-property closure join

Join DESI wedge posteriors back to FastSpecFit properties:

```bash
python workflows/sbi_inference/build_desi_wedge_property_join.py \
  --preds "${RUN_DIR}/desi_wedge_flowjax_preds.npz"
```

The join recovers `TARGETID` through `global_node_id` as an index into the source
BGS catalog, averages duplicate hemisphere copies, and writes
`desi_wedge_env_props.parquet` plus `join_report.json`. The report includes
FastSpecFit match fraction, duplicate-target count, valid-sSFR fraction, and an
RA/Dec cross-check.

Then make closure figures:

```bash
python workflows/sbi_inference/plot_property_environment_closure.py \
  --table "${RUN_DIR}/desi_wedge_env_props.parquet" \
  --outdir "${RUN_DIR}"
```

Expected qualitative closure: quenched fraction and rest-frame `g-r` increase
toward denser inferred environments, while median log sSFR decreases.

## Figure publishing

During active runs, keep run directories as the source of truth. Mirror final
figures into the canonical browsable figure root with:

```bash
scripts/sync_figures_to_canonical.sh "${RUN_DIR}" desi_wedge_flowjax_linear_si
```

The destination defaults to
`$GRAPHWEB_CANONICAL_FIGURE_DIR/desi_wedge_flowjax_linear_si`, or
`/pscratch/sd/d/dkololgi/graphweb_desi/figures/desi_wedge_flowjax_linear_si`
when the environment variable is unset.

## Common pitfalls

- Do not mix the regression wedge and path1 fiberassign wedge for
  `--abacus-gnn-arrays`; array shapes can match while scaler constants do not.
- Do not run on legacy Mpc/h DESI wedge artifacts; inference intentionally aborts
  on those metadata.
- Do not use `--scale-invariant-features` with a baseline cache/model, or a
  baseline transform with a scale-invariant cache/model.
- `SFR <= 0` is treated as quenched in closure fractions and excluded from median
  sSFR panels.
- `scripts/sync_figures_to_canonical.sh` syncs top-level image/video/HTML files
  from the source directory; nested artifacts need a separate sync command.
