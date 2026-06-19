#!/usr/bin/env python3
"""Corner plot of the inferred eigenvalue joint distributions + marginals, overlaying
Abacus truth, Abacus NPE (posterior mean) and DESI NPE (posterior mean). Shows where
in (λ1,λ2,λ3) space the cluster corner (all λ>λ_th) drains on real DESI.
"""
from __future__ import annotations
import argparse, os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

NAMES = [r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$"]


def main(a):
    apply_style()
    aself = np.load(a.abacus_self_npz)
    desi = np.load(a.desi_npz)
    rng = np.random.default_rng(0)
    dm = desi["lambda_mean"]; dm = dm[rng.choice(len(dm), min(len(dm), a.max_points), replace=False)]
    datasets = [aself["eig_truth"], aself["lambda_mean"], dm]
    labels = ["Abacus truth", "Abacus NPE", "DESI NPE"]
    colors = ["#9a9a93", ACCENT_COLORS["blue"], ACCENT_COLORS["magenta"]]
    lo, hi = a.lmin, a.lmax
    th = a.lambda_th

    fig, axes = plt.subplots(3, 3, figsize=(11, 11))
    for r in range(3):
        for c in range(3):
            ax = axes[r, c]
            if c > r:
                ax.axis("off"); continue
            if r == c:
                for X, lab, col in zip(datasets, labels, colors):
                    ax.hist(X[:, r], bins=70, range=(lo, hi), density=True, histtype="step", lw=1.8, color=col, label=lab)
                ax.axvline(th, ls=":", color="#888", lw=1)
                ax.set_yticks([])
                if r == 0:
                    ax.legend(fontsize=8, loc="upper right")
            else:
                for X, col in zip(datasets, colors):
                    x, y = X[:, c], X[:, r]
                    H, xe, ye = np.histogram2d(x, y, bins=60, range=[[lo, hi], [lo, hi]])
                    H = gaussian_filter(H, 1.3)
                    if H.max() > 0:
                        ax.contour(0.5 * (xe[:-1] + xe[1:]), 0.5 * (ye[:-1] + ye[1:]), H.T,
                                   levels=[H.max() * f for f in (0.05, 0.25, 0.6)], colors=col, linewidths=1.0)
                ax.axvline(th, ls=":", color="#888", lw=0.7); ax.axhline(th, ls=":", color="#888", lw=0.7)
                # shade the cluster corner (both axes > th)
                ax.axvspan(th, hi, ymin=(th - lo) / (hi - lo), ymax=1, color="#F5C144", alpha=0.06)
            if r == 2:
                ax.set_xlabel(NAMES[c])
            else:
                ax.set_xticklabels([])
            if c == 0 and r > 0:
                ax.set_ylabel(NAMES[r])
            elif c != 0:
                ax.set_yticklabels([])
            ax.set_xlim(lo, hi)
            if r != c:
                ax.set_ylim(lo, hi)
    frac = {lab: float(np.mean((X > th).all(1))) for X, lab in zip(datasets, labels)}
    fig.suptitle("Inferred eigenvalue joint distributions — cluster corner (gold, all λ>λ_th) "
                 f"drains on DESI\ncluster fraction: " +
                 "  ".join(f"{k}={v:.3f}" for k, v in frac.items()), fontsize=12)
    out = a.out or str(Path(a.desi_npz).parent / "eigenvalue_corner.png")
    fig.savefig(out, bbox_inches="tight", dpi=200); plt.close(fig)
    print("Saved:", out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--abacus-self-npz", required=True)
    ap.add_argument("--desi-npz", required=True)
    ap.add_argument("--out", default=None)
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--lmin", type=float, default=-1.0); ap.add_argument("--lmax", type=float, default=2.2)
    ap.add_argument("--max-points", type=int, default=20000)
    main(ap.parse_args())
