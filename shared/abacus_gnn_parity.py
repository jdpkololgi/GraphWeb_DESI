"""Abacus–DESI GNN parity helpers (coordinates, edges, training-time transforms).

Matches ``TNG/Illustris/workflows/abacus_tweb/build_abacus_sbi_cache._build_graph_from_npz``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import numpy as np
from sklearn.preprocessing import StandardScaler

CoordUnits = Literal["mpc", "mpc_per_h"]

_EDGE_LOG_EPS = 1e-6


def sky_to_xyz(
    ra_deg: np.ndarray,
    dec_deg: np.ndarray,
    z: np.ndarray,
    *,
    units: CoordUnits = "mpc",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Comoving Cartesian coords from sky position.

    Parameters
    ----------
    units
        ``mpc`` — comoving Mpc (Abacus training parity).
        ``mpc_per_h`` — legacy DESI builder (distance × h).
    """
    from astropy.cosmology import Planck18 as cosmo

    ra = np.deg2rad(np.asarray(ra_deg, dtype=np.float64))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=np.float64))
    zz = np.asarray(z, dtype=np.float64)
    dist = np.asarray(cosmo.comoving_distance(zz).value, dtype=np.float64)
    if units == "mpc_per_h":
        dist = dist * float(cosmo.h)
    x = dist * np.cos(dec) * np.cos(ra)
    y = dist * np.cos(dec) * np.sin(ra)
    zc = dist * np.sin(dec)
    return x, y, zc


def sky_to_xyz_mpc(ra_deg: np.ndarray, dec_deg: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Comoving Mpc (Abacus / training parity)."""
    return sky_to_xyz(ra_deg, dec_deg, z, units="mpc")


def duplicate_bidirectional_edges(
    edge_index: np.ndarray,
    edge_attr: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Duplicate edges with reversed direction (Abacus SBI cache convention)."""
    edge_index = np.asarray(edge_index, dtype=np.int64)
    edge_attr = np.asarray(edge_attr)
    if edge_index.shape[0] != 2:
        raise ValueError(f"edge_index must be (2, E); got {edge_index.shape}")
    rev_attr = edge_attr.copy()
    rev_attr[:, 1:4] *= -1.0
    rev_attr[:, 4] = 1.0 / np.maximum(rev_attr[:, 4], _EDGE_LOG_EPS)
    ei = np.concatenate([edge_index, edge_index[[1, 0], :]], axis=1)
    ea = np.concatenate([edge_attr, rev_attr], axis=0)
    return ei, ea


def _log_edge_cols(edge_attr: np.ndarray) -> np.ndarray:
    """Log edge length (col 0) and density contrast (col 4) in float32."""
    e = np.asarray(edge_attr, dtype=np.float32).copy()
    e[:, 0] = np.log(np.maximum(e[:, 0], _EDGE_LOG_EPS))
    e[:, 4] = np.log(np.maximum(e[:, 4], _EDGE_LOG_EPS))
    return e


def fit_edge_length_density_scaler_from_gnn_npz(
    gnn_npz: Path | str,
    *,
    make_bidirectional: bool = True,
) -> StandardScaler:
    """Fit StandardScaler on log(edge_length) and log(density_contrast) (cols 0, 4)."""
    with np.load(Path(gnn_npz).expanduser().resolve()) as data:
        edge_attr = np.asarray(data["edge_attr"], dtype=np.float32)
        if make_bidirectional:
            edge_index = np.asarray(data["edge_index"], dtype=np.int64)
            _, edge_attr = duplicate_bidirectional_edges(edge_index, edge_attr)
        else:
            edge_attr = edge_attr.copy()
    edge_attr = _log_edge_cols(edge_attr)
    scaler = StandardScaler()
    scaler.fit(edge_attr[:, [0, 4]])
    return scaler


def transform_edge_length_density(
    edge_attr: np.ndarray,
    scaler: StandardScaler,
) -> np.ndarray:
    """Apply log + StandardScaler on cols 0 and 4 (in place on a copy)."""
    e = _log_edge_cols(edge_attr)
    e[:, [0, 4]] = scaler.transform(e[:, [0, 4]]).astype(np.float32)
    return e.astype(np.float32, copy=False)


def prepare_edges_for_jraph_forward(
    edge_index: np.ndarray,
    edge_attr: np.ndarray,
    edge_scaler: StandardScaler,
    *,
    make_bidirectional: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Bidirectional duplicate (optional) then log + z-score edge cols 0, 4."""
    ei = np.asarray(edge_index, dtype=np.int64)
    ea = np.asarray(edge_attr)
    if make_bidirectional:
        ei, ea = duplicate_bidirectional_edges(ei, ea)
    ea = transform_edge_length_density(ea, edge_scaler)
    return ei, ea
