# GraphWeb_DESI Workflows and Refactor Plan

> **HISTORICAL (superseded).** This plan describes the *pre-reorg* flat layout
> and proposed a `src/graphweb_desi/` target that was **not** adopted — the repo
> moved to `workflows/` + `shared/` instead, and that migration is substantially
> complete. Kept for context only. For the current state use `ACTIVE_WORKFLOWS.md`,
> `RUNBOOK.md`, and `CLAUDE.md`.

This document maps current GraphWeb_DESI workflows and defines a practical cleanup plan that keeps existing runs reproducible while reducing script sprawl and cross-repo coupling.

## 1) Workflow map (current state)

### Workflow A: Build low-z combined catalog (**active batch script**)

**Entrypoint**

- `load_catalog.py`

**Pipeline**

1. Read many fastspec catalogs from DESI CFS.
2. Select low-z galaxy rows.
3. Join metadata + spectrophot + fastspec tables.
4. Join redshift flags and Legacy Survey photometry.
5. Write:
  - `loa-combined-lowz.fits`
  - `loa-combined-lowz-zflags.fits`
  - `loa-combined-lowz-fastspec-phot.fits`

**Notes**

- Runs at module import time (no `main()` guard).
- Contains comment/code mismatch around redshift criteria and counting message.

---

### Workflow B: Catalog filtering and coordinate transforms (**active utility class**)

**Entrypoint**

- `galaxy_catalog.py` (`GalaxyCatalog` class + local demo block)

**Pipeline**

1. Load merged FITS.
2. Apply quality and science cuts.
3. Convert sky coordinates to Cartesian positions.
4. Produce plotting/inspection outputs.

---

### Workflow C: Graph construction + feature engineering + GNN inference (**main active pipeline**)

**Entrypoint**

- `graph_catalog.py`

**Pipeline**

1. Load or rebuild cached graph objects.
2. Build alpha-complex or Delaunay graph.
3. Engineer and scale node features.
4. Load pretrained model and run inference.
5. Save VAC-like outputs and diagnostics plots.

**Dependencies**

- Local: `galaxy_catalog.py`
- Cross-repo: `TNG/Illustris/Network_stats.py` and `TNG/Illustris/Utilities.py`
- Cached artifacts in `cache/`

---

### Workflow D: Edge QA diagnostics (**active utility**)

**Entrypoint**

- `investigate_edges.py`

**Pipeline**

1. Load cached graph geometry.
2. Measure and visualize long-edge statistics.
3. Save QA figures/summary.

---

## 2) Dependency and coupling map

- `load_catalog.py` -> DESI CFS inputs
- `galaxy_catalog.py` -> merged FITS products
- `graph_catalog.py` -> `galaxy_catalog.py` + `TNG/Illustris` modules + local cache artifacts
- `investigate_edges.py` -> `cache/DESI_alpha_geom.pt`

Critical coupling:

- `graph_catalog.py` currently adds absolute `sys.path` for `TNG/Illustris`.
- Model/scaler files and cache paths are hardcoded to user-specific locations.

---

## 3) Redundancy and organization issues

1. Script/notebook duplication: logic in `graph_catalog.py` is duplicated in interactive notebooks.
2. Pipeline stages (build graph, cache, infer, plotting) are intermixed in one large script.
3. Cache naming conventions include older and newer variants, increasing ambiguity.
4. No simple repo-level runbook (`README.md`) for setup and expected outputs.

---

## 4) Refactor plan (surgical, low-risk)

## Phase 0 (1 day): lock active workflow and outputs

- Mark `graph_catalog.py` as canonical runtime pipeline.
- Declare notebooks as either `active-notebooks` or `archive-notebooks`.
- Write explicit expected outputs for each workflow.

## Phase 1 (1-2 days): centralize paths/config

- Add `config_paths.py` with defaults:
  - `DESI_CFS_ROOT`
  - `GRAPHWEB_CACHE_DIR`
  - `GRAPHWEB_OUTPUT_DIR`
  - `ILLUSTRIS_REPO_ROOT`
  - `ILLUSTRIS_MODEL_PATH`
  - `ILLUSTRIS_SCALER_PATH`
- Allow env-var overrides.
- Replace hardcoded absolute paths in:
  - `graph_catalog.py`
  - `galaxy_catalog.py`
  - `investigate_edges.py`
  - `load_catalog.py`

## Phase 2 (2-4 days): split `graph_catalog.py` by concern

- Extract modules:
  - `pipeline/graph_build.py`
  - `pipeline/features.py`
  - `pipeline/inference.py`
  - `pipeline/cache_io.py`
- Keep `graph_catalog.py` as a thin orchestrator/CLI wrapper.

## Phase 3 (1-2 days): reduce cross-repo fragility

- Replace runtime `sys.path.append(...)` with one explicit import strategy:
  - Option A: packaged dependency (`pip install -e` for Illustris tools)
  - Option B: environment variable `ILLUSTRIS_REPO_ROOT` + safe path check
- Add clear failure message when Illustris dependency is unavailable.

## Phase 4 (1-2 days): quality gates

- Add smoke tests for:
  - `load_catalog.py` argument parsing and row selection helpers
  - `graph_catalog.py` cache load path and model load path
  - `investigate_edges.py` cache existence/error handling

---

## 5) Immediate fixes before refactor

1. In `load_catalog.py`, align comments with actual cut (`0.01 <= Z <= 0.06`) or adjust logic if intended cut differs.
2. Fix counting message in `load_catalog.py` so it reports matched row count (not total scanned rows).
3. Add `if __name__ == "__main__":` guard and `main()` function to `load_catalog.py` to avoid import side effects.
4. Add CLI flags in `graph_catalog.py` for cache mode (`--use-cache/--rebuild-cache`) and graph type (`--alpha/--delaunay`).

---

## 6) Suggested target layout

- `src/graphweb_desi/`
  - `catalog/`
  - `pipeline/`
  - `inference/`
  - `config/`
- `scripts/`
  - `run_load_catalog.py`
  - `run_graph_inference.py`
  - `run_edge_qc.py`
- `notebooks/`
  - `active/`
  - `archive/`
- `cache/` (generated artifacts, optional gitignore)

Use thin scripts + reusable modules so notebook experimentation does not duplicate production workflow logic.