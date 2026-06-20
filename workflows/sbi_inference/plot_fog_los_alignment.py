#!/usr/bin/env python3
"""Proper Fingers-of-God test: recompute each galaxy's local inertia tensor (from
its Delaunay-neighbour displacements), take the MAJOR eigenvector, and measure its
alignment |cos θ| with the line of sight (radial direction). FoG elongates dense
structures along the LOS -> |cos θ| -> 1. Compares Abacus mock vs real DESI.

The earlier fog_anisotropy test used only eigenvalue magnitudes (elongation, no
direction); this adds the direction, which is what distinguishes FoG from generic
anisotropy.
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402


def xyz(ra, dec, z):
    d = cosmo.comoving_distance(np.asarray(z)).value
    r, dd = np.deg2rad(ra), np.deg2rad(dec)
    return np.vstack([d * np.cos(dd) * np.cos(r), d * np.cos(dd) * np.sin(r), d * np.sin(dd)]).T


def los_alignment(pos, edge_index):
    """|cos(major inertia axis, line of sight)| per galaxy."""
    N = len(pos)
    s, r = edge_index[0], edge_index[1]
    dvec = pos[r] - pos[s]                       # [E,3] neighbour displacement
    # symmetric inertia tensor per node (each edge contributes outer(d,d) to both ends)
    pairs = [(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]
    T = np.zeros((N, 6))
    for ci, (i, j) in enumerate(pairs):
        comp = dvec[:, i] * dvec[:, j]
        T[:, ci] = np.bincount(s, comp, N) + np.bincount(r, comp, N)
    M = np.empty((N, 3, 3))
    M[:, 0, 0], M[:, 0, 1], M[:, 0, 2] = T[:, 0], T[:, 1], T[:, 2]
    M[:, 1, 0], M[:, 1, 1], M[:, 1, 2] = T[:, 1], T[:, 3], T[:, 4]
    M[:, 2, 0], M[:, 2, 1], M[:, 2, 2] = T[:, 2], T[:, 4], T[:, 5]
    w, v = np.linalg.eigh(M)                     # ascending eigenvalues
    major = v[:, :, 2]                           # eigenvector of largest eigenvalue
    los = pos / (np.linalg.norm(pos, axis=1, keepdims=True) + 1e-9)
    return np.abs(np.sum(major * los, axis=1))   # |cos theta| in [0,1]


def main():
    apply_style()
    out = Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/fog_los_alignment.png")
    WD = "/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph"
    pre = WD + "/desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc"
    dp = np.load("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/desi_wedge_flowjax_preds.npz")
    d_edge = np.load(pre + "_gnn_arrays.npz")["edge_index"]
    d_dens = np.load(pre + "_gnn_arrays.npz")["x"][:, 2]
    d_pos = xyz(dp["ra"], dp["dec"], dp["z"])
    d_cos = los_alignment(d_pos, d_edge)

    AW = "/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3"
    at = Table.read(AW + "_wedge_targets.fits"); C = {c.upper(): c for c in at.colnames}
    a_pos = xyz(np.asarray(at[C["RA"]]), np.asarray(at[C["DEC"]]), np.asarray(at[C["Z"]]))
    a_arr = np.load(AW + "_cugraph_gnn_arrays.npz")
    a_edge = a_arr["edge_index"]; a_dens = a_arr["x"][:, 2]
    a_cos = los_alignment(a_pos, a_edge)

    # isotropy baseline: |cos| of a random axis is uniform on [0,1] (mean 0.5)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for ax, mask, ttl in [(axes[0], slice(None), "all galaxies"),
                          (axes[1], None, "densest 5% (cluster candidates)")]:
        if ttl.startswith("dense"):
            md = d_dens > np.percentile(d_dens, 95); ma = a_dens > np.percentile(a_dens, 95)
        else:
            md = np.ones(len(d_cos), bool); ma = np.ones(len(a_cos), bool)
        ax.hist(a_cos[ma], bins=40, range=(0, 1), density=True, histtype="step", lw=2, color="#9a9a93",
                label=f"Abacus (med {np.median(a_cos[ma]):.2f})")
        ax.hist(d_cos[md], bins=40, range=(0, 1), density=True, histtype="step", lw=2, color=ACCENT_COLORS["magenta"],
                label=f"DESI (med {np.median(d_cos[md]):.2f})")
        ax.axhline(1.0, color="#888", ls=":", lw=1, label="isotropic (uniform)")
        ax.set_xlabel(r"$|\cos\theta|$  (major inertia axis vs line of sight)")
        ax.set_ylabel("density"); ax.set_title(ttl); ax.legend(fontsize=9)
    fig.suptitle("Fingers-of-God test: local major-axis alignment with the line of sight")
    fig.savefig(out, bbox_inches="tight", dpi=200); plt.close(fig)
    print("Saved:", out)
    print(f"all:    Abacus median|cosθ|={np.median(a_cos):.3f}  DESI={np.median(d_cos):.3f}  (0.5=isotropic, 1=LOS-aligned)")
    md = d_dens > np.percentile(d_dens, 95); ma = a_dens > np.percentile(a_dens, 95)
    print(f"dense:  Abacus median|cosθ|={np.median(a_cos[ma]):.3f}  DESI={np.median(d_cos[md]):.3f}")
    print(f"dense frac |cosθ|>0.8 (strong LOS):  Abacus={np.mean(a_cos[ma]>0.8):.3f}  DESI={np.mean(d_cos[md]>0.8):.3f}")


if __name__ == "__main__":
    main()
