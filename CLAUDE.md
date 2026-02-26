# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository processes DESI (Dark Energy Spectroscopic Instrument) BGS (Bright Galaxy Survey) galaxy catalogs and applies graph neural network models trained on IllustrisTNG simulations to infer cosmic web environments (Void, Wall, Filament, Cluster) for observed galaxies.

## Running the Pipeline

The main inference pipeline is in `graph_catalog.py`:

```bash
python graph_catalog.py
```

This script:
1. Loads or constructs a graph from DESI BGS galaxies (Delaunay or alpha-complex)
2. Applies a pre-trained GAT model from the Illustris repository
3. Outputs predictions to `DESI_BGS_PRERELEASE_VAC.pkl`

## Architecture

### Data Flow

1. **Catalog Loading** (`load_catalog.py`)
   - Reads DESI fastspecfit catalogs from NERSC CFS: `/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/`
   - Selects galaxies with 0.01 ≤ z ≤ 0.06, SPECTYPE=GALAXY
   - Joins with redshift flags and Legacy Survey photometry
   - Outputs: `loa-combined-lowz.fits`, `loa-combined-lowz-zflags.fits`, `loa-combined-lowz-fastspec-phot.fits`

2. **Galaxy Catalog Processing** (`galaxy_catalog.py`)
   - `GalaxyCatalog` class: loads FITS files, filters by ZWARN, DELTACHI2, LOGMSTAR, BGS_TARGET
   - Converts RA/Dec to Cartesian coordinates (Mpc) using Planck18 cosmology
   - Separates galactic north/south hemispheres

3. **Graph Construction & Inference** (`graph_catalog.py`)
   - Creates `network` object (from `../TNG/Illustris/Network_stats.py`)
   - Builds alpha-complex or Delaunay graph from galaxy positions
   - Extracts node features, scales with Box-Cox transform
   - Loads pre-trained GAT model from Illustris repo
   - Outputs per-galaxy environment predictions and probabilities

### Key Dependencies

- **Illustris Repository**: `../TNG/Illustris/` provides:
  - `Network_stats.py`: Graph construction (Delaunay, alpha-complex)
  - `Utilities.py`: TNG data loading utilities (used for comparison)
  - Pre-trained model weights: `trained_gat_model_ddp_*.pth`

### Caching

Processed graphs are cached in `cache/`:
- `DESI_alpha_graph.pt` / `DESI_delaunay_graph.pt` - NetworkX graph
- `DESI_alpha_geom.pt` / `DESI_delaunay_geom.pt` - PyTorch Geometric Data object
- `DESI_alpha_features.pt` - Scaled node features
- `DESI_NETWORKalpha_zcat.pt` - Galaxy catalog with predictions

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
