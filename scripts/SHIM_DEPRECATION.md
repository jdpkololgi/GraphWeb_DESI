# Shim Deprecation Guidance

Shims are temporary root-level wrappers that preserve existing commands while
canonical scripts live under `workflows/`.

## Current shims

- `load_catalog.py`
- `graph_catalog.py`
- `galaxy_catalog.py`
- `investigate_edges.py`

## Removal criteria

Remove a shim only when all conditions are true:

1. No notebooks, docs, cron jobs, or launch scripts call the root path.
2. Users have migrated to canonical `workflows/...` commands.
3. At least one full production cycle succeeds without shim usage.
