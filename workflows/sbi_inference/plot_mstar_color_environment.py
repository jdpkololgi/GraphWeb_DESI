#!/usr/bin/env python3
"""Stellar mass vs (g-r) colour, split by inferred cosmic-web environment.

  mstar_gr_by_environment   2x2 nested KDE contours per class (reproduces the
                            reference figure: dashed enclosed-probability contours)
  mstar_gr_overlay          all four environments overlaid, one colour per class,
                            so the locus shift (redder / more massive -> denser
                            environment) is directly comparable

Uses the SI closure-join parquet (per-galaxy gr, LOGMSTAR, inferred hard_class).
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, COSMIC_WEB_COLORS  # noqa: E402

PARQUET = Path(os.environ.get("WEDGE_PARQUET",
               "/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
               "desi_wedge_flowjax_linear_si/desi_wedge_env_props.parquet"))
OUT = Path(os.environ.get("WEDGE_FIGDIR",
           "/pscratch/sd/d/dkololgi/graphweb_desi/figures/desi_wedge_flowjax_linear_si"))
CLASS_NAMES = ["Void", "Wall", "Filament", "Cluster"]   # (Sheet==Wall, Knot==Cluster)
CLASS_KEYS = ["void", "wall", "filament", "cluster"]
GX = (0.0, 2.2)     # (g-r)
GY = (9.0, 12.3)    # log M*
FRACS = [0.9, 0.7, 0.5, 0.3]   # enclosed-probability contour levels (outer->inner)
# outer->inner colour ramp that reads on the black deck theme (inner brightest)
RAMP = ["#2f6f5e", "#2aa198", "#6b8cff", "#FFFFFF"]


def load():
    df = pd.read_parquet(PARQUET)
    df = df[df["has_props"]].copy()
    df = df[np.isfinite(df["gr"]) & np.isfinite(df["LOGMSTAR"])]
    df = df[(df["gr"] > GX[0]) & (df["gr"] < GX[1]) &
            (df["LOGMSTAR"] > GY[0]) & (df["LOGMSTAR"] < GY[1])]
    return df


def kde_grid(x, y, n=140):
    xx, yy = np.mgrid[GX[0]:GX[1]:n*1j, GY[0]:GY[1]:n*1j]
    kde = gaussian_kde(np.vstack([x, y]))
    Z = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    return xx, yy, Z


def enclosed_levels(Z, fracs):
    """Density levels enclosing the given probability fractions."""
    flat = np.sort(Z.ravel())[::-1]
    cs = np.cumsum(flat); cs /= cs[-1]
    lvls = [flat[min(np.searchsorted(cs, f), len(flat) - 1)] for f in fracs]
    return lvls


def _axfmt(ax):
    ax.set_xlim(GX); ax.set_ylim(GY)
    ax.grid(True, ls="--", alpha=0.15)


def fig_by_env(df):
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 9.0), sharex=True, sharey=True,
                             layout="constrained")
    for k, ax in enumerate(axes.ravel()):
        sub = df[df["hard_class"] == k]
        xx, yy, Z = kde_grid(sub["gr"].to_numpy(), sub["LOGMSTAR"].to_numpy())
        lvls = enclosed_levels(Z, FRACS)   # FRACS 0.9->0.3 gives increasing density
        # keep (level, colour) pairs strictly increasing for contour()
        pairs = sorted(zip(lvls, RAMP), key=lambda t: t[0])
        lv = [p[0] for p in pairs]; cols = [p[1] for p in pairs]
        keep = [0] + [i for i in range(1, len(lv)) if lv[i] > lv[i-1]]
        ax.contour(xx, yy, Z, levels=[lv[i] for i in keep],
                   colors=[cols[i] for i in keep],
                   linestyles="--", linewidths=1.8)
        ax.text(0.96, 0.94, f"{CLASS_NAMES[k]}", transform=ax.transAxes,
                ha="right", va="top", fontsize=17, color=COSMIC_WEB_COLORS[CLASS_KEYS[k]])
        ax.text(0.96, 0.85, f"N={len(sub):,}", transform=ax.transAxes,
                ha="right", va="top", fontsize=11, color="#B8B8B8")
        _axfmt(ax)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$(g-r)$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\log_{10}(M_*/M_\odot)$")
    _save(fig, "mstar_gr_by_environment")


def fig_overlay(df):
    fig, ax = plt.subplots(figsize=(7.6, 6.6), layout="constrained")
    for k in range(4):
        sub = df[df["hard_class"] == k]
        xx, yy, Z = kde_grid(sub["gr"].to_numpy(), sub["LOGMSTAR"].to_numpy())
        lvls = sorted(enclosed_levels(Z, [0.9, 0.5]))   # 2 nested contours
        c = COSMIC_WEB_COLORS[CLASS_KEYS[k]]
        ax.contour(xx, yy, Z, levels=lvls, colors=[c, c],
                   linestyles=["--", "-"], linewidths=[1.4, 2.2], alpha=0.95)
        ax.plot([], [], color=c, lw=2.4, label=f"{CLASS_NAMES[k]}  (N={len(sub):,})")
    ax.legend(fontsize=12, loc="lower right", frameon=False)
    ax.set_xlabel(r"$(g-r)$"); ax.set_ylabel(r"$\log_{10}(M_*/M_\odot)$")
    ax.set_title("Mass–colour by inferred environment (50% & 90% contours)", fontsize=14)
    _axfmt(ax)
    _save(fig, "mstar_gr_overlay")


def fig_all_colormass(df):
    """Whole-wedge colour-mass diagram + marginal (g-r) histogram — the figure
    where the classic blue-cloud / red-sequence bimodality actually shows for
    this sample (SFR-M* is single-peaked because FastSpecFit sSFR is)."""
    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(7.8, 7.4))
    gs = GridSpec(2, 1, height_ratios=[1, 4], hspace=0.04, figure=fig)
    axh = fig.add_subplot(gs[0]); ax = fig.add_subplot(gs[1], sharex=axh)
    # main colour-mass hexbin
    hb = ax.hexbin(df["gr"], df["LOGMSTAR"], gridsize=60, extent=[*GX, *GY],
                   mincnt=1, cmap="viridis", bins="log")
    ax.set_xlim(GX); ax.set_ylim(GY)
    ax.set_xlabel(r"$(g-r)$"); ax.set_ylabel(r"$\log_{10}(M_*/M_\odot)$")
    ax.grid(True, ls="--", alpha=0.15)
    # marginal colour histogram (the bimodality: blue cloud + red sequence)
    axh.hist(df["gr"], bins=np.arange(GX[0], GX[1], 0.03), color="#8ab6ff",
             histtype="stepfilled", alpha=0.85)
    axh.axvline(0.80, ls=":", color="#F2F2F2", alpha=0.6)  # blue/red divide
    axh.text(0.30, 0.7, "blue\ncloud", transform=axh.transAxes, ha="center",
             va="center", color="#8ab6ff", fontsize=11)
    axh.text(0.62, 0.7, "red\nsequence", transform=axh.transAxes, ha="center",
             va="center", color="#e07a7a", fontsize=11)
    axh.set_ylabel("N", fontsize=12); axh.tick_params(labelbottom=False)
    axh.set_title("Mass–colour bimodality (whole wedge, N={:,})".format(len(df)),
                  fontsize=13)
    cb = fig.colorbar(hb, ax=[axh, ax], fraction=0.05, pad=0.02); cb.set_label(r"$N_{\rm gal}$")
    _save(fig, "mstar_gr_all_bimodal")


def _save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", OUT / f"{stem}.*", flush=True)


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    df = load()
    print(f"N (valid gr & M*): {len(df):,}", flush=True)
    for k in range(4):
        print(f"  {CLASS_NAMES[k]:9s} N={int((df['hard_class']==k).sum()):,}"
              f"  median M*={df.loc[df['hard_class']==k,'LOGMSTAR'].median():.2f}"
              f"  median (g-r)={df.loc[df['hard_class']==k,'gr'].median():.3f}", flush=True)
    fig_by_env(df)
    fig_overlay(df)
    fig_all_colormass(df)


if __name__ == "__main__":
    main()
