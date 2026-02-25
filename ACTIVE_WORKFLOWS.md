# Active Workflow Index

This file is the Phase 0 quick reference for what to run in this repository now.
For pscratch organization and migration env vars, see `/global/homes/d/dkololgi/PSCRATCH_LAYOUT.md`.

## Active

- Catalog assembly:
  - Canonical: `workflows/catalog/load_catalog.py`
  - Compatibility shim: `load_catalog.py`
- Graph inference pipeline:
  - Canonical: `workflows/graph_inference/graph_catalog.py`
  - Compatibility shim: `graph_catalog.py`
- Utility:
  - Canonical: `workflows/utilities/galaxy_catalog.py`
  - Canonical: `workflows/utilities/investigate_edges.py`
  - Compatibility shims: `galaxy_catalog.py`, `investigate_edges.py`

## Experimental

- Notebook-driven variants under this repo are exploratory unless explicitly promoted.

## Legacy/To retire

- Any duplicate notebook/script flow superseded by `graph_catalog.py`.

## Compatibility policy

- Root-level script names are temporary wrappers to keep existing commands working.
- Prefer calling canonical workflow paths for new scripts, docs, and SLURM launchers.
