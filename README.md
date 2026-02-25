# GraphWeb_DESI

GraphWeb_DESI builds a graph from DESI BGS galaxies and applies a pretrained GAT model (from `TNG/Illustris`) to infer cosmic web environments.

## Main entrypoint

- Canonical: `workflows/graph_inference/graph_catalog.py`
- Compatibility wrapper: `graph_catalog.py`

This script now has an explicit CLI and no longer runs large ad-hoc plotting blocks by default.

## Quick start

Run from repo root:

```bash
python graph_catalog.py
```

Canonical path form:

```bash
python workflows/graph_inference/graph_catalog.py
```

Other entrypoints follow the same pattern:
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

Defaults come from `config_paths.py`. Useful env vars include:

- `GRAPHWEB_CACHE_DIR`
- `GRAPHWEB_VAC_OUTPUT_PATH`
- `ILLUSTRIS_REPO_ROOT`
- `ILLUSTRIS_GAT_MODEL_PATH`
- `ILLUSTRIS_SCALER_PATH`
- `TNG_REFERENCE_CATALOG_PATH`

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
