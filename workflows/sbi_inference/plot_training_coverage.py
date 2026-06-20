#!/usr/bin/env python3
"""Training-coverage test for the cluster deficit. The GNN was trained on Abacus
node features (box-cox scaled, scaler fit on the Abacus train split). DESI features
get the SAME scaler. If DESI's dense / cluster-candidate galaxies land beyond the
Abacus training support in the scaled features, the model is extrapolating — which
would explain why it under-predicts their collapse strength (clusters -> filaments).

Quantifies, per node feature: the systematic DESI shift vs training, and the
fraction of DESI (and DESI cluster-candidate) galaxies above the Abacus training
p99.9 / max. Figure overlays the scaled distributions on a log-y axis.
"""
from __future__ import annotations
import os, pickle, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

NAMES = ["Degree", "Clustering", "Density", "NeighDensity", "I_eig1", "I_eig2", "I_eig3"]
SHOW = [0, 2, 3]  # density-related features (clusters live here)


def main():
    apply_style()
    CACHE = "/pscratch/sd/d/dkololgi/abacus/sbi_caches/path1_flowjax_3d_lineareig/processed_jraph_data_mc1e+09_v2_scaled_3_linear_eig.pkl"
    with open(CACHE, "rb") as f:
        cache = pickle.load(f)
    nfs = cache["node_feature_scaler"]
    ab = np.asarray(cache["graph"].nodes)               # Abacus, box-cox scaled (all nodes)
    train = np.asarray(cache["masks"][0])
    ab_tr = ab[train]                                   # training support

    WD = "/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph"
    pre = WD + "/desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc"
    desi_raw = np.load(pre + "_gnn_arrays.npz")["x"].astype(np.float64)
    desi = nfs.transform(desi_raw + 1e-6)               # same box-cox as inference
    hard = np.load("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/desi_wedge_flowjax_preds.npz")["hard_class"]
    clu = hard == 3

    print("feature        train_p99.9   DESI_p99.9   DESI_max   frac_DESI>train_p99.9   frac_DESI>train_max   (cluster-cand frac>p99.9)")
    dens = desi_raw[:, 2]; cand = dens > np.percentile(dens, 95)   # dense cluster candidates
    for k in range(7):
        tp = np.percentile(ab_tr[:, k], 99.9); tmax = ab_tr[:, k].max()
        f99 = np.mean(desi[:, k] > tp); fmax = np.mean(desi[:, k] > tmax)
        fcand = np.mean(desi[cand, k] > tp)
        flag = "  <== OOD" if (f99 > 0.02 or fmax > 0.005) else ""
        print(f"{NAMES[k]:13s}  {tp:8.2f}    {np.percentile(desi[:,k],99.9):8.2f}   {desi[:,k].max():7.2f}   {f99:8.4f}            {fmax:8.4f}        {fcand:8.4f}{flag}")

    fig, axes = plt.subplots(1, len(SHOW), figsize=(6 * len(SHOW), 5))
    for ax, k in zip(axes, SHOW):
        lo = min(ab_tr[:, k].min(), desi[:, k].min()); hi = max(ab_tr[:, k].max(), desi[:, k].max())
        bins = np.linspace(lo, hi, 80)
        ax.hist(ab_tr[:, k], bins=bins, density=True, histtype="step", lw=2, color="#9a9a93", label="Abacus train")
        ax.hist(desi[:, k], bins=bins, density=True, histtype="step", lw=2, color=ACCENT_COLORS["blue"], label="DESI (all)")
        ax.hist(desi[clu, k], bins=bins, density=True, histtype="step", lw=2, color=ACCENT_COLORS["magenta"], label="DESI inferred cluster")
        ax.axvline(ab_tr[:, k].max(), color="#F5C144", ls="--", lw=1.5, label="Abacus train max")
        ax.set_yscale("log"); ax.set_xlabel(f"{NAMES[k]} (box-cox scaled)"); ax.set_ylabel("density")
        ax.set_title(NAMES[k])
        if k == SHOW[0]:
            ax.legend(fontsize=8)
    fig.suptitle("Training-coverage check: DESI dense-feature tail vs Abacus training support")
    out = "/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/training_coverage.png"
    fig.savefig(out, bbox_inches="tight", dpi=200); plt.close(fig)
    print("Saved:", out)


if __name__ == "__main__":
    main()
