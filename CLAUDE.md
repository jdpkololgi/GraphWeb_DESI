# CLAUDE.md

Guidance for Claude Code / automation agents working in this repository.
(`AGENTS.md` is a symlink to this file.)

## Start here

- `README.md` — repository orientation.
- `ACTIVE_WORKFLOWS.md` — current canonical entrypoints and their status.
- `RUNBOOK.md` — validated launch commands and environment setup.
- `~/.claude/CLAUDE.md` — cross-repo map, conda envs, and Perlmutter job recipes.
- the `nersc` skill — NERSC/Slurm depth (auto-loads for Perlmutter work).

## Project overview

Processes DESI BGS (Bright Galaxy Survey) galaxy catalogs and applies graph
neural network models trained on IllustrisTNG (in `../TNG/Illustris/`) to infer
cosmic web environments (Void / Wall / Filament / Cluster) for observed galaxies.

## Repository layout (post-reorg)

Canonical code lives under `workflows/` and `shared/`. Root-level scripts
(`graph_catalog.py`, `load_catalog.py`, `galaxy_catalog.py`,
`investigate_edges.py`, `config_paths.py`) are **thin compatibility shims** —
prefer the canonical paths (see `scripts/SHIM_DEPRECATION.md`).

| Path | Purpose | Env |
| --- | --- | --- |
| `workflows/catalog/load_catalog.py` | Assemble low-z combined DESI catalogs from CFS fastspecfit. | `cosmic_env` + `desienv` |
| `workflows/utilities/galaxy_catalog.py` | `GalaxyCatalog`: quality/mass cuts, RA/Dec/Z → Cartesian (Planck18), N/S split. | `cosmic_env` |
| `workflows/graph_construction/build_desi_bgs_gudhi_graph.py` | Build DESI BGS Gudhi alpha/Delaunay graph artifacts (CPU, memory-heavy). | `cosmic_env` |
| `workflows/graph_construction/desi_graph_features_cugraph.py` | GPU/cuGraph node-feature extraction. | `rapids-gnn`, GPU node |
| `workflows/graph_construction/subset_desi_graph_wedge.py` | Subset a graph to an RA/Dec/z wedge. | `cosmic_env` |
| `workflows/graph_inference/graph_catalog.py` | **Main inference driver:** graph → features → GAT → VAC. | `cosmic_env` |
| `workflows/utilities/investigate_edges.py` | Long-edge QA diagnostics. | `cosmic_env` |
| `workflows/visualization/` | 3D wedge cosmic-web notebooks (+ archived ones). | `cosmic_env` |
| `shared/config_paths.py` | Env-var-driven path configuration. | — |

## Running the pipeline

See `RUNBOOK.md` for exact commands. Canonical inference entrypoint:

```bash
python workflows/graph_inference/graph_catalog.py --help
```

It loads or builds a DESI graph (alpha-complex or Delaunay), engineers and
Box-Cox-scales node features, applies the pre-trained GAT, and writes per-galaxy
predictions to `DESI_BGS_PRERELEASE_VAC.pkl`. Processed graphs are cached under
`cache/` (`DESI_alpha_*`, `DESI_delaunay_*`).

cuGraph feature extraction is the GPU-heavy step — see the Perlmutter recipes in
`~/.claude/CLAUDE.md` and the `.cursor/rules/rapids-gnn-env.mdc` rule.

## Cross-repo dependency

`graph_catalog.py` imports `Network_stats` / `Utilities` from
`../TNG/Illustris/` and loads its pre-trained GAT weights
(`trained_gat_model_ddp_*.pth`). The graph machinery is shared so DESI and TNG
point clouds run through identical construction/feature code.

## Environment classification

Four T-Web classes: `0 = Void`, `1 = Wall`, `2 = Filament`, `3 = Cluster`.

## Data locations on NERSC

- Fastspecfit catalogs: `/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/`
- Redshift catalog: `/global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits`
- Legacy Survey photometry: `/global/cfs/cdirs/desi/vac/dr2/lsdr9-photometry/loa/v1.0/observed-targets`
