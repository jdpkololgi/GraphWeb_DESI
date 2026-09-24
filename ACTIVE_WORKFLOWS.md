# Active Workflow Index

## P12-A Loa application (2026-09-24)

Active plan: `../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md`.
Implement the frozen U-PATCH/P12-A manifest, Loa input/response crosswalk and
golden-mock replay before the bounded real-data trial. No validated new inference
entrypoint is asserted yet. Legacy GAT/Gudhi/Jraph workflows below remain usable
for their historical scope; they do not implement this new posterior handoff.
Slurm, reserved phase-access and science-release gates remain in force.

This file is the Phase 0 quick reference for what to run in this repository now.
For scratch layout, conda envs, and Perlmutter job recipes, see `~/.claude/CLAUDE.md`;
for launch commands see `RUNBOOK.md`.

## Active

- Low-z catalog assembly for the GAT classifier:
  - Canonical: `workflows/catalog/load_catalog.py`
  - Compatibility shim: `load_catalog.py`
- Bright BGS catalog assembly for Gudhi/Jraph:
  - Canonical: `workflows/catalog/build_bgs_maglim_catalog.py`
  - Produces the magnitude-limited FITS table used by the full-graph Gudhi stack.
- GAT graph inference pipeline:
  - Canonical: `workflows/graph_inference/graph_catalog.py`
  - Compatibility shim: `graph_catalog.py`
- Gudhi/cuGraph graph construction:
  - Full graph: `workflows/graph_construction/build_desi_bgs_gudhi_graph.py`
  - Feature export: `workflows/graph_construction/desi_graph_features_cugraph.py`
  - Wedge subset: `workflows/graph_construction/subset_desi_graph_wedge.py`
  - Use Mpc coordinates (`--coord-units mpc`, default) for Abacus/Jraph parity.
- Jraph inference:
  - Active DESI wedge inference: `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`
  - Legacy PyG-cache wedge builder: `workflows/jraph_inference/build_desi_wedge_jraph_cache.py`
  - Feature-parity experiment: `workflows/jraph_inference/experiment_desi_feature_parity.py`
- Utility:
  - Canonical: `workflows/utilities/galaxy_catalog.py`
  - Canonical: `workflows/utilities/investigate_edges.py`
  - Compatibility shims: `galaxy_catalog.py`, `investigate_edges.py`
- Shared helpers:
  - Config paths: `shared/config_paths.py` (shim: `config_paths.py`)
  - Abacus/Jraph coordinate and edge-feature parity: `shared/abacus_gnn_parity.py`

## Visualization

- Active: `workflows/visualization/visualize_desi_wedge_cweb_3d.ipynb`
- Archived notebooks: `workflows/visualization/archive/` (LOA catalog, graph subvolume; see README)

## Experimental

- Notebook-driven variants under this repo are exploratory unless explicitly promoted.

## Legacy/To retire

- Any duplicate notebook/script flow superseded by `graph_catalog.py`.
- Older Mpc/h DESI wedge artifacts are deprecated; prefer the expanded Mpc-parity
  paths in `workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt`.

## Compatibility policy

- Root-level script names are temporary wrappers to keep existing commands working.
- Prefer calling canonical workflow paths for new scripts, docs, and SLURM launchers.
- Config paths canonical module: `shared/config_paths.py` (shim: `config_paths.py`)
