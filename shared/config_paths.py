"""Centralized path configuration for GraphWeb_DESI workflows.

Phase 1.5 goal:
- Preserve current defaults.
- Add canonical pscratch layout targets for gradual migration.
"""

from __future__ import annotations

import os
from pathlib import Path


def _env(name: str, default: str) -> str:
    return os.environ.get(name, default)


GRAPHWEB_REPO_ROOT = _env(
    "GRAPHWEB_REPO_ROOT",
    str(Path(__file__).resolve().parents[1]),
)
DK_SCRATCH_ROOT = _env("DK_SCRATCH_ROOT", "/pscratch/sd/d/dkololgi")
GRAPHWEB_SCRATCH_ROOT = _env("GRAPHWEB_SCRATCH_ROOT", f"{DK_SCRATCH_ROOT}/graphweb_desi")
GRAPHWEB_CACHE_DIR = _env(
    "GRAPHWEB_CACHE_DIR",
    f"{GRAPHWEB_REPO_ROOT}/cache",
)
GRAPHWEB_OUTPUT_DIR = _env(
    "GRAPHWEB_OUTPUT_DIR",
    GRAPHWEB_REPO_ROOT,
)
GRAPHWEB_VAC_OUTPUT_PATH = _env(
    "GRAPHWEB_VAC_OUTPUT_PATH",
    f"{GRAPHWEB_OUTPUT_DIR}/DESI_BGS_PRERELEASE_VAC.pkl",
)
# Catalogs live on pscratch (moved off the home quota 2026-07-15; home hit 100% and broke writes).
# NOTE: the old default pointed at {GRAPHWEB_REPO_ROOT}/loa-combined-lowz.fits, which was already
# stale (the files had been under data/). Canonical home is now GRAPHWEB_CATALOG_DIR on pscratch.
GRAPHWEB_CATALOG_DIR = _env(
    "GRAPHWEB_CATALOG_DIR",
    f"{GRAPHWEB_SCRATCH_ROOT}/catalogs",
)
GRAPHWEB_CATALOG_PATH = _env(
    "GRAPHWEB_CATALOG_PATH",
    f"{GRAPHWEB_CATALOG_DIR}/loa-combined-lowz.fits",
)
GRAPHWEB_CATALOG_ZFLAGS_PATH = _env(
    "GRAPHWEB_CATALOG_ZFLAGS_PATH",
    f"{GRAPHWEB_CATALOG_DIR}/loa-combined-lowz-zflags.fits",
)
GRAPHWEB_CATALOG_FASTSPEC_PATH = _env(
    "GRAPHWEB_CATALOG_FASTSPEC_PATH",
    f"{GRAPHWEB_CATALOG_DIR}/loa-combined-lowz-fastspec-phot.fits",
)

# Canonical pscratch layout (opt-in via env vars in current migration stage)
GRAPHWEB_CANONICAL_CACHE_DIR = _env(
    "GRAPHWEB_CANONICAL_CACHE_DIR",
    f"{GRAPHWEB_SCRATCH_ROOT}/cache",
)
GRAPHWEB_CANONICAL_OUTPUT_DIR = _env(
    "GRAPHWEB_CANONICAL_OUTPUT_DIR",
    f"{GRAPHWEB_SCRATCH_ROOT}/outputs",
)
GRAPHWEB_CANONICAL_FIGURE_DIR = _env(
    "GRAPHWEB_CANONICAL_FIGURE_DIR",
    f"{GRAPHWEB_SCRATCH_ROOT}/figures",
)

ILLUSTRIS_REPO_ROOT = _env(
    "ILLUSTRIS_REPO_ROOT",
    "/global/homes/d/dkololgi/TNG/Illustris",
)
ILLUSTRIS_SCALER_PATH = _env(
    "ILLUSTRIS_SCALER_PATH",
    f"{ILLUSTRIS_REPO_ROOT}/features_scaler.pkl",
)
ILLUSTRIS_GAT_MODEL_PATH = _env(
    "ILLUSTRIS_GAT_MODEL_PATH",
    f"{ILLUSTRIS_REPO_ROOT}/trained_gat_model_ddp_2026-01-15.pth",
)

TNG_REFERENCE_CATALOG_PATH = _env(
    "TNG_REFERENCE_CATALOG_PATH",
    "/pscratch/sd/d/dkololgi/TNG300-1",
)

DESI_FASTSPEC_CATALOGS_DIR = _env(
    "DESI_FASTSPEC_CATALOGS_DIR",
    "/global/cfs/cdirs/desi/vac/dr2/fastspecfit/loa/v1.0/catalogs/",
)
DESI_ZCAT_FILE = _env(
    "DESI_ZCAT_FILE",
    "/global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits",
)
DESI_TRACTORPHOT_DIR = _env(
    "DESI_TRACTORPHOT_DIR",
    "/global/cfs/cdirs/desi/vac/dr2/lsdr9-photometry/loa/v1.0/observed-targets",
)
