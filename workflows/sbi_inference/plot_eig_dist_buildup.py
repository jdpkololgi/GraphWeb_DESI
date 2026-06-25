#!/usr/bin/env python3
"""Build-up version of the 3-way eigenvalue-distribution figure for the talk.

Emits two PNGs with SHARED bins + axis limits so they overlay cleanly across
successive slides:
  slide 1  eig_dist_buildup_abacus_only.png  -> Abacus truth + Abacus NPE
  slide 2  eig_dist_buildup_with_desi.png    -> + DESI NPE posterior means

Same theme/colours as plot_desi_wedge_flowjax.py's eigenvalue_distributions_3way:
truth = grey fill, Abacus NPE = blue step, DESI NPE = magenta step.
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

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

_GREY = "#9a9a93"


def _draw(out_path, abacus_truth, abacus_npe, desi_npe, lam_th, bins, xlim, title):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for k in range(3):
        ax = axes[k]
        ax.hist(abacus_truth[:, k], bins=bins[k], density=True, alpha=0.4,
                color=_GREY, label="Abacus truth")
        ax.hist(abacus_npe[:, k], bins=bins[k], density=True, histtype="step", lw=2,
                color=ACCENT_COLORS["blue"], label="Abacus NPE")
        if desi_npe is not None:
            ax.hist(desi_npe[:, k], bins=bins[k], density=True, histtype="step", lw=2,
                    color=ACCENT_COLORS["magenta"], label="DESI NPE")
        ax.axvline(lam_th, ls=":", color="#F2F2F2", alpha=0.6)
        ax.set_xlim(xlim[k])
        ax.set_xlabel(rf"$\lambda_{k+1}$"); ax.set_ylabel("density" if k == 0 else "")
        if k == 0:
            ax.legend(fontsize=9)
    fig.suptitle(title)
    fig.savefig(out_path, bbox_inches="tight", dpi=200); plt.close(fig)
    print("Saved:", out_path, flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abacus-self-npz", type=Path, required=True)
    ap.add_argument("--desi-preds-npz", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--lambda-threshold", type=float, default=0.2)
    args = ap.parse_args()
    apply_style()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    aself = np.load(args.abacus_self_npz)
    abacus_truth = aself["eig_truth"]; abacus_npe = aself["lambda_mean"]
    desi_npe = np.load(args.desi_preds_npz)["lambda_mean"]

    # shared bins + xlim across both figures (driven by all three series)
    bins, xlim = [], []
    for k in range(3):
        allv = np.concatenate([abacus_truth[:, k], abacus_npe[:, k], desi_npe[:, k]])
        lo, hi = np.percentile(allv, 0.5), np.percentile(allv, 99.5)
        bins.append(np.linspace(lo, hi, 80)); xlim.append((lo, hi))

    _draw(args.out_dir / "eig_dist_buildup_abacus_only.png",
          abacus_truth, abacus_npe, None, args.lambda_threshold, bins, xlim,
          "Eigenvalue distributions: Abacus truth vs Abacus NPE (posterior mean)")
    _draw(args.out_dir / "eig_dist_buildup_with_desi.png",
          abacus_truth, abacus_npe, desi_npe, args.lambda_threshold, bins, xlim,
          "Eigenvalue distributions: Abacus truth vs Abacus NPE vs DESI NPE (posterior mean)")


if __name__ == "__main__":
    main()
