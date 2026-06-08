# GraphWeb_DESI

GraphWeb_DESI contains DESI BGS graph workflows for cosmic-web inference. The
original production path builds a DESI graph and applies a pretrained PyTorch GAT
model from `TNG/Illustris` to infer four environment classes. Newer wedge work
uses Gudhi/cuGraph graph features and an Abacus-trained Jraph regression model to
predict T-Web eigenvalues on DESI sky cuts.

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
