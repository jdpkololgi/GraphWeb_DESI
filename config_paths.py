"""Centralized path configuration for GraphWeb_DESI workflows.

Defaults preserve current behavior and can be overridden by environment vars.
"""

from __future__ import annotations

import os
from pathlib import Path


def _env(name: str, default: str) -> str:
    return os.environ.get(name, default)


GRAPHWEB_REPO_ROOT = _env(
    "GRAPHWEB_REPO_ROOT",
    str(Path(__file__).resolve().parent),
)
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
GRAPHWEB_CATALOG_PATH = _env(
    "GRAPHWEB_CATALOG_PATH",
    f"{GRAPHWEB_REPO_ROOT}/loa-combined-lowz.fits",
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
