#!/usr/bin/env python3
"""Vertical (3-row) eigenvalue-distribution figure for the DR3-KP talk slide 29.

A tall, narrow variant of plot_eig_dist_buildup.py's 3-way figure, meant to sit
in a LEFT column next to the pipeline diagram on the right of the slide. Panels
stack lambda1 (top) / lambda2 / lambda3 (bottom).

FONT SIZE / KEYNOTE 35 pt
-------------------------
matplotlib font sizes are points (1 pt = 1/72 inch) and are ABSOLUTE in the
saved file. We fix figsize and save WITHOUT bbox_inches="tight" so the output
PDF's physical size == figsize exactly. Placed in Keynote at that native size
(no rescale), fontsize=35 renders as 35 pt. If you scale the figure to height H
inches on the slide, effective pt = 35 * (H / FIG_H). FIG_H is printed at the end.

truth = grey fill, Abacus NPE = blue step, DESI NPE = magenta step (same theme).
"""
from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

_GREY = "#9a9a93"

# Physical size of the figure, in inches, = the size to place it at in Keynote.
# Sized for a left column on a 16:9 slide (13.333 x 7.5 in): ~half width, most height.
FIG_W, FIG_H = 7.0, 7.3
FS = 25        # Keynote pt for axis labels + tick numbers
FS_LEG = 22    # legend a touch smaller so 3 entries fit the width without colliding


def _draw(out_stem, abacus_truth, abacus_npe, desi_npe, lam_th, bins, xlim):
    fig, axes = plt.subplots(3, 1, figsize=(FIG_W, FIG_H), layout="constrained")
    handles = None
    for k in range(3):
        ax = axes[k]
        h_truth = ax.hist(abacus_truth[:, k], bins=bins[k], density=True, alpha=0.4,
                          color=_GREY, label="Abacus truth")
        ax.hist(abacus_npe[:, k], bins=bins[k], density=True, histtype="step", lw=2.5,
                color=ACCENT_COLORS["blue"], label="Abacus NPE")
        if desi_npe is not None:
            ax.hist(desi_npe[:, k], bins=bins[k], density=True, histtype="step", lw=2.5,
                    color=ACCENT_COLORS["magenta"], label="DESI NPE")
        ax.axvline(lam_th, ls=":", color="#F2F2F2", alpha=0.6, lw=2)
        ax.set_xlim(xlim[k])
        ax.set_xlabel(rf"$\lambda_{k+1}$", fontsize=FS, labelpad=1)
        ax.xaxis.set_major_locator(MaxNLocator(3))
        ax.yaxis.set_major_locator(MaxNLocator(3))
        for lbl in ax.get_xticklabels() + ax.get_yticklabels():
            lbl.set_fontsize(FS)
        ax.set_ylabel("")
        if k == 0:
            handles, labels = ax.get_legend_handles_labels()
    # single legend as a horizontal strip along the top, out of the data
    fig.legend(handles, labels, fontsize=FS_LEG, loc="outside upper center",
               ncol=3, frameon=False, handlelength=1.1, columnspacing=1.0,
               handletextpad=0.4)
    fig.supylabel("density", fontsize=FS)
    for ext in ("pdf", "png"):
        fig.savefig(f"{out_stem}.{ext}", dpi=300)  # NO bbox_inches -> size == figsize
    plt.close(fig)
    print(f"Saved: {out_stem}.pdf / .png  (native {FIG_W}x{FIG_H} in, labels/ticks {FS} pt)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abacus-self-npz", type=Path,
                    default=Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
                                 "abacus_self_linear_si/abacus_self_flowjax_preds.npz"))
    ap.add_argument("--desi-preds-npz", type=Path,
                    default=Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
                                 "desi_wedge_flowjax_linear_si/desi_wedge_flowjax_preds.npz"))
    ap.add_argument("--out-dir", type=Path,
                    default=Path("/pscratch/sd/d/dkololgi/graphweb_desi/figures/"
                                 "desi_wedge_flowjax_linear_si"))
    ap.add_argument("--lambda-threshold", type=float, default=0.2)
    args = ap.parse_args()
    apply_style()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    aself = np.load(args.abacus_self_npz)
    abacus_truth = aself["eig_truth"]; abacus_npe = aself["lambda_mean"]
    desi_npe = np.load(args.desi_preds_npz)["lambda_mean"]

    bins, xlim = [], []
    for k in range(3):
        allv = np.concatenate([abacus_truth[:, k], abacus_npe[:, k], desi_npe[:, k]])
        lo, hi = np.percentile(allv, 0.5), np.percentile(allv, 99.5)
        bins.append(np.linspace(lo, hi, 80)); xlim.append((lo, hi))

    _draw(args.out_dir / "eig_dist_vertical_abacus_only",
          abacus_truth, abacus_npe, None, args.lambda_threshold, bins, xlim)
    _draw(args.out_dir / "eig_dist_vertical_3way",
          abacus_truth, abacus_npe, desi_npe, args.lambda_threshold, bins, xlim)


if __name__ == "__main__":
    main()
