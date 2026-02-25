# GraphWeb Workflow Reorganization

This note records the workflow-folder migration for `GraphWeb_DESI`.

## Canonical workflow paths

- `workflows/catalog/load_catalog.py`
- `workflows/graph_inference/graph_catalog.py`
- `workflows/utilities/galaxy_catalog.py`
- `workflows/utilities/investigate_edges.py`

## Compatibility shims

Root-level scripts are kept as wrappers so existing commands continue to work:

- `load_catalog.py`
- `graph_catalog.py`
- `galaxy_catalog.py`
- `investigate_edges.py`

Each shim forwards to the corresponding canonical module using `runpy`.

## Import resolution

Canonical scripts add repo-root path resolution so imports from `config_paths.py`
and other root-level modules still resolve when scripts are launched from nested
workflow directories.

## Usage policy

- Prefer canonical `workflows/...` paths in new docs and launch scripts.
- Keep shims during migration; remove only after no active callers depend on them.
