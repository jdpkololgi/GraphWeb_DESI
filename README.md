# GraphWeb_DESI

GraphWeb_DESI contains DESI BGS graph workflows for cosmic-web inference. The
original production path builds a DESI graph and applies a pretrained PyTorch GAT
model from `TNG/Illustris` to infer four environment classes. Newer wedge work
uses Gudhi/cuGraph graph features with Abacus-trained models: a Jraph regression
model for point estimates of T-Web eigenvalues, and a FlowJAX/SBI neural
posterior estimator for per-galaxy eigenvalue posteriors on DESI sky cuts.

## Main entrypoint

- Canonical: `workflows/graph_inference/graph_catalog.py`
- Compatibility wrapper: `graph_catalog.py`

This script now has an explicit CLI and no longer runs large ad-hoc plotting blocks by default.

## Pipeline Map

### A. GAT Environment Classification

Use `workflows/graph_inference/graph_catalog.py` for the all-sky DESI BGS VAC-style
classification workflow. It loads or constructs an alpha-complex/Delaunay graph,
extracts node features via the Illustris graph utilities, applies the pretrained
GAT checkpoint, and writes class probabilities for:

- 0: void
- 1: wall
- 2: filament
- 3: cluster

### B. Jraph Eigenvalue Regression

Use this stack for Abacus-parity DESI wedge inference rather than VAC production:

1. Build a bright BGS FITS catalog with `workflows/catalog/build_bgs_maglim_catalog.py`.
2. Build a hemisphere-split Gudhi graph with `workflows/graph_construction/build_desi_bgs_gudhi_graph.py`.
3. Export Abacus-style node and edge features with `workflows/graph_construction/desi_graph_features_cugraph.py`.
4. Induce an RA/Dec/z wedge with `workflows/graph_construction/subset_desi_graph_wedge.py`.
5. Run DESI inference with `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`.

Details and canonical Perlmutter paths live in
`workflows/catalog/BUILD_JRAPH_CACHES.md` and
`workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt`.

The validated Jraph path uses comoving Mpc coordinates (`--coord-units mpc`,
the graph-builder default) for Abacus training parity. Older Mpc/h wedge
artifacts are deprecated; see `workflows/catalog/to-delete_README.md`.

### C. FlowJAX/SBI Posterior Inference

Use `workflows/sbi_inference/infer_desi_wedge_flowjax.py` when the science
question needs uncertainty-aware T-Web posteriors rather than a deterministic
regression label. This path consumes the same Mpc-parity DESI wedge graph arrays
as the Jraph stack, loads the Abacus-trained FlowJAX/GNN model from the sibling
Illustris repo (`ILLUSTRIS_ROOT`), and writes per-galaxy posterior means,
posterior widths, class probabilities, and hard classes.

Primary companion scripts:

- `workflows/sbi_inference/plot_desi_wedge_flowjax.py` for truth-free DESI
  diagnostics and comparison plots.
- `workflows/sbi_inference/infer_abacus_self_flowjax.py` for Abacus self
  inference references used in DESI-vs-Abacus overlays.
- `workflows/sbi_inference/build_desi_wedge_property_join.py` and
  `plot_property_environment_closure.py` for joining inferred environments to
  LOA FastSpecFit galaxy properties.
- `scripts/sync_figures_to_canonical.sh` for mirroring finished figures into the
  canonical figure tree.

See `workflows/sbi_inference/README.md` for the validated command sequence,
model/cache constraints, and common failure modes.

## Quick start

Run from repo root:

```bash
python graph_catalog.py
```

Canonical path form:

```bash
python workflows/graph_inference/graph_catalog.py
```

Other GAT-stack entrypoints follow the same pattern:
- `workflows/catalog/load_catalog.py` (shim: `load_catalog.py`)
- `workflows/utilities/galaxy_catalog.py` (shim: `galaxy_catalog.py`)
- `workflows/utilities/investigate_edges.py` (shim: `investigate_edges.py`)

Default behavior:
- graph type: `alpha`
- cache mode: `rebuild`
- writes VAC to `GRAPHWEB_VAC_OUTPUT_PATH` (from config)
- shows a summary histogram plot

## Common run modes

Use cache if available, otherwise build:

```bash
python graph_catalog.py --cache-mode prefer-cache
```

Fail if cache files are missing:

```bash
python graph_catalog.py --cache-mode cache-only
```

Use Delaunay graph instead of alpha-complex:

```bash
python graph_catalog.py --graph-type delaunay
```

Run without summary plot:

```bash
python graph_catalog.py --no-summary-plot
```

## CLI flags

Key options:
- `--graph-type {alpha,delaunay}`
- `--cache-mode {rebuild,prefer-cache,cache-only}`
- `--first-moment-matching`
- `--cache-dir <path>`
- `--model-path <path>`
- `--scaler-path <path>`
- `--vac-output-path <path>`
- `--reference-catalog-path <path>`
- `--no-summary-plot`

## Path and environment configuration

Defaults come from `shared/config_paths.py` (shim: `config_paths.py`). Useful
env vars for the GAT classification stack include:

- `GRAPHWEB_CACHE_DIR`
- `GRAPHWEB_VAC_OUTPUT_PATH`
- `ILLUSTRIS_REPO_ROOT`
- `ILLUSTRIS_GAT_MODEL_PATH`
- `ILLUSTRIS_SCALER_PATH`
- `TNG_REFERENCE_CATALOG_PATH`

Jraph inference loads model code from the Illustris training repo with
`ILLUSTRIS_ROOT` instead. In the current NERSC layout this usually points to the
same directory as `ILLUSTRIS_REPO_ROOT`.

FlowJAX/SBI inference also uses `ILLUSTRIS_ROOT`; it must resolve before this
repo's local `shared/` package so that the trained GNN encoder, FlowJAX helpers,
and eigenvalue transforms match the Abacus training code.

Example:

```bash
export GRAPHWEB_CACHE_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/cache
export GRAPHWEB_VAC_OUTPUT_PATH=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/DESI_BGS_PRERELEASE_VAC.pkl
python graph_catalog.py --cache-mode prefer-cache
```

## Outputs

Primary output:
- VAC pickle with predicted environment class and class probabilities.

Cached artifacts (in `GRAPHWEB_CACHE_DIR`):
- graph object (`*_graph.pt`)
- torch geometric data (`*_geom.pt`)
- scaled features (`*_features.pt`)
- zcat snapshot (`*_zcat.pt`)

SBI/FlowJAX wedge outputs live under a run directory, usually
`$GRAPHWEB_SCRATCH_ROOT/flowjax_inference_outputs/<run-name>/`:
- `desi_wedge_flowjax_preds.npz` with `lambda_mean`, `lambda_std`,
  `classprob`, `hard_class`, `p_exceed`, sky coordinates, and `global_node_id`.
- `summary.json` with provenance, scaler constants, class fractions, and
  ordering/consistency diagnostics.
- Figure products from `plot_desi_wedge_flowjax.py` and the closure-test scripts.
