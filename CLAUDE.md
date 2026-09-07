# CLAUDE.md

Guidance for Claude Code / automation agents working in this repository.
(`AGENTS.md` is a symlink to this file.)

## Start here

- `~/TNG/Illustris/SCIENCE_LOG.md` — **read first**: current scientific direction, open threads, and recent decisions shared between Desktop and NERSC agents. Lives in the Illustris repo; pull before reading.
- `CONTEXT.md` — durable GraphWeb_DESI context. If it conflicts with the shared Illustris science log, follow the science log.
- `README.md` — repository orientation.
- `ACTIVE_WORKFLOWS.md` — current canonical entrypoints and their status.
- `RUNBOOK.md` — validated launch commands and environment setup.
- `~/.claude/CLAUDE.md` — cross-repo map, conda envs, and Perlmutter job recipes.
- the `nersc` skill — NERSC/Slurm depth (auto-loads for Perlmutter work).

## Project Overview

This repository processes DESI (Dark Energy Spectroscopic Instrument) BGS
(Bright Galaxy Survey) galaxy catalogs and applies graph neural network models
trained on simulations to infer cosmic-web environments for observed galaxies.
It currently contains two analysis paths: a PyTorch GAT classifier for
VAC-style environment labels and a Gudhi/cuGraph/Jraph wedge workflow for
Abacus-parity eigenvalue regression.

## Running the Pipeline

The main GAT classification pipeline is in
`workflows/graph_inference/graph_catalog.py`:

```bash
python workflows/graph_inference/graph_catalog.py
```

This script:
1. Loads or constructs a graph from DESI BGS galaxies (Delaunay or alpha-complex)
2. Applies a pre-trained GAT model from the Illustris repository
3. Outputs predictions to `DESI_BGS_PRERELEASE_VAC.pkl`

The root-level `graph_catalog.py` is a compatibility shim. Prefer canonical
`workflows/...` paths in new commands and docs.

The Jraph wedge workflow is a separate path:

1. `workflows/catalog/build_bgs_maglim_catalog.py`
2. `workflows/graph_construction/build_desi_bgs_gudhi_graph.py`
3. `workflows/graph_construction/desi_graph_features_cugraph.py`
4. `workflows/graph_construction/subset_desi_graph_wedge.py`
5. `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`

The FlowJAX/SBI path reuses steps 1–4, then
`workflows/sbi_inference/infer_desi_wedge_flowjax.py` (see
`workflows/sbi_inference/README.md` and `RUNBOOK.md`).

Use `--coord-units mpc` (the default) in the Gudhi graph builder for Abacus
training parity. The cuGraph feature step runs in the `rapids-gnn` env on a GPU
node; everything else uses `cosmic_env` (see `~/.claude/CLAUDE.md` and
`.cursor/rules/rapids-gnn-env.mdc`). On Slurm, call the env’s absolute python
after unsetting `PYTHONPATH` (`.cursor/rules/conda-env-srun-python-path.mdc`).

## Architecture

### Data Flow

1. **Low-z Catalog Loading** (`workflows/catalog/load_catalog.py`)
   - Reads DESI fastspecfit catalogs from NERSC CFS: `/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/`
   - Selects galaxies with 0.01 ≤ z ≤ 0.06, SPECTYPE=GALAXY from a **hard-coded**
     FastSpecFit filename list (not a directory glob)
   - Joins with redshift flags (`ZCAT_PRIMARY`) and Legacy Survey photometry
   - Outputs: `loa-combined-lowz.fits`, `loa-combined-lowz-zflags.fits`, `loa-combined-lowz-fastspec-phot.fits`
   - **Location (moved 2026-07-15):** these live on pscratch at
     `/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/` — NOT in the repo. Home hit its 40 GiB quota
     (100%) and broke writes/commits, so ~1 GB of catalogs was moved off. Resolve paths via
     `shared/config_paths.py`: `GRAPHWEB_CATALOG_DIR`, `GRAPHWEB_CATALOG_PATH`,
     `GRAPHWEB_CATALOG_ZFLAGS_PATH`, `GRAPHWEB_CATALOG_FASTSPEC_PATH` (all env-overridable).
     Do not re-add large data under the repo/home.
   - Note: this low-z product covers **z ≈ 0.01–0.06 only** — it does NOT overlap the GraphWeb-BGS
     VAC range (0.15–0.55), so it cannot supply M_abs/ABSMAG for that work.
   - `--zcat-file` defaults to `DESI_ZCAT_FILE`
     (`/global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits`).

2. **Bright maglim catalog** (`workflows/catalog/build_bgs_maglim_catalog.py`)
   - Streaming zall filter: ZWARN==0, DELTACHI2>=25, SPECTYPE=GALAXY, BGS_BRIGHT bits
   - No redshift or stellar-mass cut; columns are `TARGET_RA`/`TARGET_DEC` (not `RA`/`DEC`)
   - Default `--zall-path` is the **public** DR2 zall under `desi/public/dr2/...`, not `DESI_ZCAT_FILE`

3. **Galaxy Catalog Processing** (`workflows/utilities/galaxy_catalog.py`)
   - `GalaxyCatalog` class: loads FITS files, filters by ZWARN, DELTACHI2, LOGMSTAR>=9, BGS_TARGET
   - Converts `RA`/`Dec` to Cartesian coordinates (Mpc) using Planck18 cosmology
   - Separates galactic north/south hemispheres and reorders the table to `[north, south]`
   - GAT-stack utility for `loa-combined-lowz.fits`. Do not feed it the maglim catalog.
   - GAT inference itself uses Illustris `network(masscut=9.0, from_DESI=True)` /
     `DESI_GAL_CAT`, then `SimpleGAT` with 10 node features.

4. **GAT Graph Construction & Inference** (`workflows/graph_inference/graph_catalog.py`)
   - Creates a DESI-aware Illustris `network` object
   - Builds alpha-complex or Delaunay graph from galaxy positions
   - Extracts node features, scales with Box-Cox transform (or `--first-moment-matching`)
   - Loads pre-trained GAT model from Illustris repo
   - Outputs per-galaxy environment predictions and probabilities

5. **Jraph Wedge Inference** (`workflows/jraph_inference/`)
   - Starts from a bright BGS zall-derived catalog with no redshift or stellar-mass cuts
   - Builds a hemisphere-split Gudhi graph in comoving Mpc (**no north–south edges**)
   - Exports seven node features and five edge features with Abacus-style names
   - Subsets an RA/Dec/z wedge (catalog row order must match the parent graph)
   - Outputs predicted λ1, λ2, λ3 eigenvalues, T-Web-like classes, and diagnostic plots

6. **FlowJAX/SBI Wedge Posteriors** (`workflows/sbi_inference/`)
   - Same Mpc-parity DESI wedge GNN arrays as Jraph, plus path1 Abacus edge scaler
   - Per-galaxy posterior samples over ordered eigenvalues; classes from `lambda_th=0.2`
   - See `workflows/sbi_inference/README.md`

### Key Dependencies

- **Illustris Repository**: `/global/homes/d/dkololgi/TNG/Illustris` by default
  - GAT path uses `ILLUSTRIS_REPO_ROOT` from `shared/config_paths.py`
  - Jraph path uses `ILLUSTRIS_ROOT` to import `shared/graph_net_models.py`
  - FlowJAX/SBI also uses `ILLUSTRIS_ROOT` (put it on `sys.path` before this repo)
  - Pre-trained GAT model weights: `trained_gat_model_ddp_*.pth`
  - Abacus-trained Jraph checkpoints live under the Illustris/Abacus run outputs
- **Path config**: use `shared/config_paths.py`; root-level `config_paths.py` is a shim.

### Caching

GAT processed graphs are cached in `cache/` by default:
- `DESI_alpha_graph.pt` / `DESI_delaunay_graph.pt` - NetworkX graph
- `DESI_alpha_geom.pt` / `DESI_delaunay_geom.pt` - PyTorch Geometric Data object
- `DESI_alpha_features.pt` - Scaled node features
- `DESI_NETWORKalpha_zcat.pt` - Galaxy catalog with predictions

Jraph/Gudhi artifacts are typically written under
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/` and documented in
`workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt`.

## Environment Classification

Four cosmic web classes (from T-Web formalism):
- 0: Void
- 1: Wall
- 2: Filament
- 3: Cluster

## Data Locations on NERSC

- Fastspecfit catalogs: `/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/`
- Redshift catalog: `/global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits`
- Legacy Survey photometry: `/global/cfs/cdirs/desi/vac/dr2/lsdr9-photometry/loa/v1.0/observed-targets`

## graphify

This project has a knowledge graph at graphify-out/ with god nodes, community structure, and cross-file relationships.

**Sibling repo (cross-project):** `~/TNG/Illustris` — simulation-side ML pipelines (Abacus T-Web, Jraph regression, SBI, GCN). This repo consumes its outputs: `shared/graph_net_models.py`, GAT model weights, and Abacus-trained Jraph checkpoints via `ILLUSTRIS_ROOT` / `ILLUSTRIS_REPO_ROOT`.

Rules:
- MANDATORY: Before using Read, Grep, Glob, or Bash to explore the codebase, run graphify first.
- For codebase questions, run `graphify query "<question>"`. Use `graphify path "<A>" "<B>"` for relationships and `graphify explain "<concept>"` for focused concepts. These return a scoped subgraph, usually much smaller than GRAPH_REPORT.md or raw grep output.
- If graphify-out/wiki/index.md exists, use it for broad navigation instead of raw source browsing.
- Read graphify-out/GRAPH_REPORT.md only for broad architecture review or when query/path/explain do not surface enough context.
- After modifying code, run `graphify update .` to keep the graph current (AST-only, no API cost); then refresh the global with `graphify global add graphify-out/graph.json --as GraphWeb_DESI`.
- Only use Read/Grep/Glob directly when graphify has already oriented you and you need to modify or debug specific lines.
- This rule applies to subagents too — include it in every subagent prompt involving code exploration.

**Cross-repo questions** (GraphWeb_DESI ↔ Illustris, model imports, training parity, shared paths): use the global graph:
- `graphify query "<question>" --graph ~/.graphify/global-graph.json`
- `graphify path "<A>" "<B>" --graph ~/.graphify/global-graph.json`
- `graphify explain "<concept>" --graph ~/.graphify/global-graph.json`
