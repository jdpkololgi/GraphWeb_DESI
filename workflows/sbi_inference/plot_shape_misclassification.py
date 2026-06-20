#!/usr/bin/env python3
"""Proof test for the shape/morphology hypothesis of the cluster deficit.

Among DENSE galaxies (cluster candidates), does the inferred-cluster rate fall with
local anisotropy (elongation)? And is the rate-vs-anisotropy curve the SAME for
Abacus and DESI? If yes, the cluster deficit is explained by DESI simply having
MORE elongated dense structures (a distribution shift the model handles
consistently). If DESI sits lower at fixed anisotropy, there is an extra domain
shift beyond shape. An Oaxaca-style decomposition splits the two.
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402


def aniso(x):
    e = np.sort(np.abs(x[:, 4:7]) + 1e-9, 1); return e[:, 2] / e[:, 0]


def main():
    apply_style()
    out = "/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/shape_misclassification.png"
    WD = "/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph"
    dx = np.load(WD + "/desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_gnn_arrays.npz")["x"]
    dh = np.load("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/desi_wedge_flowjax_preds.npz")["hard_class"]
    asf = np.load("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear/abacus_self_flowjax_preds.npz")
    axf = np.load("/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz")["x"]
    ax = axf[asf["node_index"]]; ah = asf["hard_class"]

    # dense = top 10% by Density, per dataset
    md = dx[:, 2] > np.percentile(dx[:, 2], 90); ma = ax[:, 2] > np.percentile(ax[:, 2], 90)
    da, aa = aniso(dx)[md], aniso(ax)[ma]
    dc, ac = (dh[md] == 3).astype(float), (ah[ma] == 3).astype(float)
    print(f"dense galaxies: Abacus n={ma.sum()} cluster-rate={ac.mean():.3f} | DESI n={md.sum()} cluster-rate={dc.mean():.3f}")

    # shared anisotropy bins (quantiles of the pooled dense anisotropy)
    edges = np.quantile(np.concatenate([aa, da]), np.linspace(0, 1, 11))
    edges[-1] += 1e-6
    cen = 0.5 * (edges[:-1] + edges[1:])
    def curve(av, cv):
        idx = np.digitize(av, edges) - 1; idx = np.clip(idx, 0, len(cen) - 1)
        r = np.array([cv[idx == b].mean() if np.any(idx == b) else np.nan for b in range(len(cen))])
        p = np.array([np.mean(idx == b) for b in range(len(cen))])
        return r, p
    rA, pA = curve(aa, ac); rD, pD = curve(da, dc)

    # Oaxaca decomposition of the dense cluster-rate gap (Abacus - DESI)
    valid = ~np.isnan(rA) & ~np.isnan(rD)
    rateA = np.nansum(pA * np.where(np.isnan(rA), 0, rA))
    rateD = np.nansum(pD * np.where(np.isnan(rD), 0, rD))
    # counterfactual: DESI anisotropy distribution x Abacus curve
    cf = np.nansum(pD[valid] * rA[valid]) / pD[valid].sum()
    gap = rateA - rateD
    shape_part = rateA - cf          # due to DESI being more elongated (distribution)
    resid_part = cf - rateD          # due to the curve differing at fixed shape
    print(f"dense cluster-rate gap (Abacus-DESI) = {gap:.3f}")
    print(f"  explained by DESI being MORE elongated (distribution shift): {shape_part:.3f} ({100*shape_part/gap:.0f}%)")
    print(f"  residual (model treats DESI differently at fixed shape):     {resid_part:.3f} ({100*resid_part/gap:.0f}%)")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    axes[0].hist(aa, bins=np.linspace(1, 15, 50), density=True, histtype="step", lw=2, color="#9a9a93", label=f"Abacus dense (med {np.median(aa):.1f})")
    axes[0].hist(da, bins=np.linspace(1, 15, 50), density=True, histtype="step", lw=2, color=ACCENT_COLORS["magenta"], label=f"DESI dense (med {np.median(da):.1f})")
    axes[0].set_xlabel("anisotropy $I_{eig,max}/I_{eig,min}$"); axes[0].set_ylabel("density")
    axes[0].set_title("Dense-galaxy elongation"); axes[0].legend(fontsize=9)
    axes[1].plot(cen, rA, "-o", color="#9a9a93", ms=4, label="Abacus")
    axes[1].plot(cen, rD, "-o", color=ACCENT_COLORS["magenta"], ms=4, label="DESI")
    axes[1].set_xlabel("anisotropy (elongation)"); axes[1].set_ylabel("P(inferred cluster | dense)")
    axes[1].set_title("Cluster rate vs shape — same curve = deficit is just shape")
    axes[1].legend(fontsize=9)
    fig.suptitle("Is the cluster deficit explained by dense-structure shape? "
                 f"(shape {100*shape_part/gap:.0f}% / residual {100*resid_part/gap:.0f}%)")
    fig.savefig(out, bbox_inches="tight", dpi=200); plt.close(fig)
    print("Saved:", out)


if __name__ == "__main__":
    main()
