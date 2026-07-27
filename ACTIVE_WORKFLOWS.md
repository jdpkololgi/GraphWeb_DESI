# Active Workflow Index

This file is the Phase 0 quick reference for what to run in this repository now.
For scratch layout, conda envs, and Perlmutter job recipes, see `~/.claude/CLAUDE.md`;
for launch commands see `RUNBOOK.md`.

## Active

- Low-z catalog assembly for the GAT classifier:
  - Canonical: `workflows/catalog/load_catalog.py`
  - Compatibility shim: `load_catalog.py`
  - Writes to `GRAPHWEB_CATALOG_DIR` on pscratch by default (not the repo/home).
- Bright BGS catalog assembly for Gudhi/Jraph/SBI:
  - Canonical: `workflows/catalog/build_bgs_maglim_catalog.py`
  - Produces the magnitude-limited FITS table used by the full-graph Gudhi stack.
- GAT graph inference pipeline:
  - Canonical: `workflows/graph_inference/graph_catalog.py`
  - Compatibility shim: `graph_catalog.py`
- Gudhi/cuGraph graph construction:
  - Full graph: `workflows/graph_construction/build_desi_bgs_gudhi_graph.py`
  - Feature export: `workflows/graph_construction/desi_graph_features_cugraph.py`
  - Wedge subset: `workflows/graph_construction/subset_desi_graph_wedge.py`
  - Use Mpc coordinates (`--coord-units mpc`, default) for Abacus/Jraph/SBI parity.
- Jraph inference:
  - Active DESI wedge inference: `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`
  - Legacy PyG-cache wedge builder: `workflows/jraph_inference/build_desi_wedge_jraph_cache.py`
  - Feature-parity experiment: `workflows/jraph_inference/experiment_desi_feature_parity.py`
- FlowJAX/SBI posterior inference:
  - DESI wedge NPE inference: `workflows/sbi_inference/infer_desi_wedge_flowjax.py`
  - Abacus self-inference reference: `workflows/sbi_inference/infer_abacus_self_flowjax.py`
  - DESI posterior diagnostics: `workflows/sbi_inference/plot_desi_wedge_flowjax.py`
  - Property-closure join/plots:
    `workflows/sbi_inference/build_desi_wedge_property_join.py`,
    `workflows/sbi_inference/plot_property_environment_closure.py`
  - Workflow README: `workflows/sbi_inference/README.md`
  - Requires the same Mpc-parity DESI wedge products, path1 Abacus edge-scaler
    arrays, and `ILLUSTRIS_ROOT` for the Illustris FlowJAX model code.
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
- Candidate luminosity-channel support diagnostic:
  - `workflows/sbi_inference/desi_absmag_kcorr.py`
  - Computes DESI `ABSMAG_RP1` with the DESI LSS k+e implementation and prints
    a marginal comparison with a raw Abacus cut-sky sample.
  - This is not a production pass/fail gate: the current DESI and Abacus
    selections, redshift distributions, and magnitude conventions are not
    matched. See `RUNBOOK.md` before interpreting its output.
- Additional SBI investigation scripts under `workflows/sbi_inference/`
  (coverage, FoG, velocity dispersion, smoothing-scale, mass-anchored recovery,
  skewer animations, etc.) are research diagnostics, not canonical launch paths.

## Legacy/To retire

- Any duplicate notebook/script flow superseded by `graph_catalog.py`.
- Older Mpc/h DESI wedge artifacts are deprecated; prefer the expanded Mpc-parity
  paths in `workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt`.

## Compatibility policy

- Root-level script names are temporary wrappers to keep existing commands working.
- Prefer calling canonical workflow paths for new scripts, docs, and SLURM launchers.
- Config paths canonical module: `shared/config_paths.py` (shim: `config_paths.py`)
