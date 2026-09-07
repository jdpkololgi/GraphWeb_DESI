# GraphWeb Workflow Reorganization

This note records the workflow-folder migration for `GraphWeb_DESI`.

## Canonical paths

- `workflows/catalog/load_catalog.py`
- `workflows/graph_inference/graph_catalog.py`
- `workflows/utilities/galaxy_catalog.py`
- `workflows/utilities/investigate_edges.py`
- `shared/config_paths.py`

## Compatibility shims

Root-level scripts are kept as wrappers so existing commands continue to work:

- `load_catalog.py`
- `graph_catalog.py`
- `galaxy_catalog.py`
- `investigate_edges.py`
- `config_paths.py`

Script shims forward to the corresponding canonical module using `runpy`.
`config_paths.py` re-exports `shared.config_paths` for old imports.

## Import resolution

Canonical scripts add repo-root path resolution so imports from shared modules
still resolve when scripts are launched from nested workflow directories.
Use `shared/config_paths.py` for new code; root-level `config_paths.py` is a
temporary shim.

## Post-migration subsystems

Newer workflow groups were added after the original four-script migration:

- `workflows/graph_construction/` for Gudhi graph building, cuGraph features,
  and wedge subsetting.
- `workflows/jraph_inference/` for Abacus-trained Jraph DESI wedge inference and
  experiments.
- `workflows/sbi_inference/` for FlowJAX NPE posteriors, property joins, and
  transfer diagnostics (README there is the detailed companion).
- `workflows/visualization/` for the path1 Jraph 3D notebook (plus `archive/`).
- `shared/` for config, Abacus/Jraph parity helpers, and plot style.
- `scripts/` / `templates/job.sbatch` for NERSC helpers and a minimal sbatch
  template.

Catalog note: `load_catalog.py` (low-z FastSpecFit, `RA`/`DEC`) and
`build_bgs_maglim_catalog.py` (bright zall, `TARGET_RA`/`TARGET_DEC`) are
different products. Root shims exist only for the original four scripts plus
`config_paths.py`.

## Usage policy

- Prefer canonical `workflows/...` paths in new docs and launch scripts.
- Keep shims during migration; remove only after no active callers depend on them.
