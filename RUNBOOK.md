# GraphWeb_DESI Runbook

This runbook lists validated workflow entrypoints and launch commands.

## Environment setup

Recommended shell setup for GraphWeb workflows:

```bash
source ~/.bashrc
conda activate cosmic_env
desienv
```

If `desimodel` import remains unavailable in your shell context, use the
canonical script paths and environment used in validation notes. The catalog
workflow also includes a hardcoded `desimodel` fallback path.

## Canonical workflows

### Catalog assembly

Canonical:

```bash
python workflows/catalog/load_catalog.py --help
```

Compatibility wrapper:

```bash
python load_catalog.py --help
```

Example run:

```bash
python workflows/catalog/load_catalog.py --out-dir .
```

### Graph inference pipeline

Canonical:

```bash
python workflows/graph_inference/graph_catalog.py --help
```

Compatibility wrapper:

```bash
python graph_catalog.py --help
```

Example no-plot cache-only probe:

```bash
python workflows/graph_inference/graph_catalog.py --cache-mode cache-only --no-summary-plot
```

If model checkpoint path differs from config default:

```bash
python workflows/graph_inference/graph_catalog.py --model-path /path/to/trained_gat_model.pth
```

### Gudhi/cuGraph/Jraph wedge inference

This path is for DESI wedge eigenvalue regression with an Abacus-trained Jraph
model. It is not a replacement command for the all-sky GAT VAC workflow above.
Run it on Perlmutter from an environment with DESI, Gudhi/cuGraph where needed,
JAX/Haiku/Jraph, and access to the Illustris training repo.

Build the bright BGS catalog if it does not already exist:

```bash
python workflows/catalog/build_bgs_maglim_catalog.py \
  --out-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits
```

Build the full Mpc-parity graph and export graph features:

```bash
python workflows/graph_construction/build_desi_bgs_gudhi_graph.py \
  --catalog-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc \
  --alpha-sq inf \
  --coord-units mpc

python workflows/graph_construction/desi_graph_features_cugraph.py \
  --metadata-path /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_delaunay_metadata.json \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc \
  --out-prefix desi_bgs_cugraph
```

Subset the expanded Abacus-parity DESI wedge and run inference:

```bash
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

python workflows/graph_construction/subset_desi_graph_wedge.py \
  --graph-metadata /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_delaunay_metadata.json \
  --catalog-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --parent-gnn-arrays /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_bgs_cugraph_gnn_arrays.npz \
  --parent-gnn-metadata /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_bgs_cugraph_gnn_metadata.json \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph \
  --out-prefix desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc \
  --ra-min 120 --ra-max 160 --dec-min 14.5 --dec-max 30.6 --z-min 0.2 --z-max 0.3

python workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py \
  --abacus-run-dir "${ABACUS_RUN_DIR}" \
  --calibration-cache "${CALIBRATION_CACHE}" \
  --abacus-gnn-arrays "${ABACUS_WEDGE_GNN_ARRAYS}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --output-dir "$(dirname "${INFER_DIR}")" \
  --run-name "${INFER_RUN_NAME}"
```

`INFER_DIR` in `workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt` is the full
intended run directory, so the command passes its parent as `--output-dir` and
the basename as `--run-name`. Jraph inference uses `ILLUSTRIS_ROOT` to find
`shared/graph_net_models.py`; the GAT classifier uses `ILLUSTRIS_REPO_ROOT` from
`shared/config_paths.py`.

### FlowJAX/SBI DESI wedge posterior inference

This path runs an Abacus-trained FlowJAX neural posterior estimator on the DESI
wedge. It returns per-galaxy posterior eigenvalue summaries and class
probabilities, not just a hard class or point-regression output. Use it after the
Mpc-parity Gudhi/cuGraph/Jraph wedge products above exist.

Detailed notes and caveats:

```bash
less workflows/sbi_inference/README.md
```

Validated Perlmutter launcher:

```bash
bash workflows/sbi_inference/run_infer_desi_wedge_flowjax.sh
```

Manual command shape:

```bash
export ILLUSTRIS_ROOT=/global/homes/d/dkololgi/TNG/Illustris
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

FLOWJAX_RUN=/pscratch/sd/d/dkololgi/abacus/sbi_runs/path1_wedge_flowjax_3d_Bcorrected_linear
MODEL=${FLOWJAX_RUN}/flowjax_sbi_model_seed_42_<timestamp>.pkl
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
  --run-name desi_wedge_flowjax_linear_si \
  --scale-invariant-features
```

`--scale-invariant-features` must match a scale-invariant training cache and
model. For transfer diagnostics only, the script also exposes
`--edge-domain-adapt` and `--node-domain-adapt`.

Generate DESI diagnostics after inference:

```bash
RUN_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear_si

python workflows/sbi_inference/plot_desi_wedge_flowjax.py \
  --preds-npz "${RUN_DIR}/desi_wedge_flowjax_preds.npz" \
  --summary-json "${RUN_DIR}/summary.json" \
  --output-dir "${RUN_DIR}" \
  --no-regression
```

Join FastSpecFit properties and make closure figures:

```bash
python workflows/sbi_inference/build_desi_wedge_property_join.py \
  --preds "${RUN_DIR}/desi_wedge_flowjax_preds.npz"

python workflows/sbi_inference/plot_property_environment_closure.py \
  --table "${RUN_DIR}/desi_wedge_env_props.parquet" \
  --outdir "${RUN_DIR}"
```

Mirror final top-level figures into the canonical figure root:

```bash
scripts/sync_figures_to_canonical.sh "${RUN_DIR}" desi_wedge_flowjax_linear_si
```

### Utility workflows

Canonical:

```bash
python workflows/utilities/galaxy_catalog.py --help
python workflows/utilities/investigate_edges.py --help
```

Compatibility wrappers:

```bash
python galaxy_catalog.py --help
python investigate_edges.py --help
```

## Notes

- Active workflows and statuses: `ACTIVE_WORKFLOWS.md`
- Reorg notes: `docs/migration/WORKFLOW_REORG.md`
- Shim policy: `scripts/SHIM_DEPRECATION.md`
- Prefer canonical `workflows/...` paths in new docs and launch scripts.
