#!/usr/bin/env python3
"""Star-forming main-sequence (SFR vs M*) diagrams for the DESI wedge, split by
inferred cosmic-web environment, plus a continuous version colour-coded by the
posterior-mean tidal eigenvalues.

Figures (PNG + PDF, dark deck theme):
  sfms_all            whole wedge, hexbin counts, MS + Green-Valley/Red-Sequence
  sfms_by_environment 2x2 hexbin counts, one panel per inferred class
  sfms_eigen_continuous  1x3 hexbin coloured by mean lambda_1 / _2 / _3 per cell

Regions: main sequence fit to the star-forming population (log sSFR > -10.7);
Blue Cloud above MS-0.6 dex, Green Valley MS-0.6..MS-1.2, Red Sequence below.
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, COSMIC_WEB_COLORS  # noqa: E402

_CLASS_KEYS = ["void", "wall", "filament", "cluster"]

PARQUET = Path(os.environ.get("WEDGE_PARQUET",
               "/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
               "desi_wedge_flowjax_linear_si/desi_wedge_env_props.parquet"))
OUT = Path(os.environ.get("WEDGE_FIGDIR",
           "/pscratch/sd/d/dkololgi/graphweb_desi/figures/desi_wedge_flowjax_linear_si"))
CLASS_NAMES = ["Void", "Wall", "Filament", "Cluster"]

# plotting window
MX = (8.5, 12.6)   # log M*
MY = (-3.4, 2.2)   # log SFR
GV_UP, GV_LO = 0.8, 1.6   # dex below MS ridge: green-valley band edges (brackets
                          # the CIGALE bimodal trough / top of the red sequence)


def load():
    df = pd.read_parquet(PARQUET)
    df = df[df["has_ssfr"] & np.isfinite(df["LOGMSTAR"]) & np.isfinite(df["SFR"])].copy()
    df["logSFR"] = np.log10(df["SFR"].to_numpy())
    # drop FastSpecFit SFR-floor artefacts (log sSFR spikes far below the locus)
    df = df[(df["log_sSFR"] > -13.5) & np.isfinite(df["logSFR"])]
    return df


def fit_ms(df):
    # Fit the star-forming RIDGE: binned median log SFR of clearly-star-forming
    # galaxies (a plain polyfit over all SF galaxies flattens because the blue
    # cloud thins/bends at high mass and picks up green-valley contamination).
    sf = df[df["log_sSFR"] > -10.5]
    x = sf["LOGMSTAR"].to_numpy(); y = sf["logSFR"].to_numpy()
    edges = np.arange(9.4, 11.4, 0.2)
    xc, ym = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        mm = (x >= lo) & (x < hi)
        if mm.sum() >= 30:
            xc.append(0.5 * (lo + hi)); ym.append(np.median(y[mm]))
    m, b = np.polyfit(xc, ym, 1)
    return m, b


def draw_regions(ax, m, b, shade=True):
    x = np.array(MX)
    ms = m * x + b
    if shade:
        ax.fill_between(x, ms - GV_UP, MY[1], color="#4f7bd0", alpha=0.07, zorder=0)
        ax.fill_between(x, ms - GV_LO, ms - GV_UP, color="#3fae6b", alpha=0.10, zorder=0)
        ax.fill_between(x, MY[0], ms - GV_LO, color="#d05050", alpha=0.07, zorder=0)
    ax.plot(x, ms, "--", color="#F2F2F2", lw=1.8, zorder=5, label="Main Sequence")
    ax.plot(x, ms - GV_UP, ":", color="#3fae6b", lw=1.6, zorder=5)
    ax.plot(x, ms - GV_LO, ":", color="#d05050", lw=1.6, zorder=5)


def _axfmt(ax):
    ax.set_xlim(MX); ax.set_ylim(MY)
    ax.set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
    ax.set_ylabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")


def fig_all(df, m, b):
    fig, ax = plt.subplots(figsize=(7.2, 5.6), layout="constrained")
    draw_regions(ax, m, b)
    hb = ax.hexbin(df["LOGMSTAR"], df["logSFR"], gridsize=55, extent=[*MX, *MY],
                   mincnt=1, cmap="viridis", bins="log")
    cb = fig.colorbar(hb, ax=ax); cb.set_label(r"$N_{\rm gal}$")
    ax.text(0.04, 0.93, "Blue Cloud", transform=ax.transAxes, color="#7ea6ff", fontsize=12, va="top")
    ax.text(0.30, 0.42, "Green Valley", transform=ax.transAxes, color="#57c281", fontsize=11,
            rotation=np.degrees(np.arctan(m)), va="center")
    ax.text(0.72, 0.10, "Red Sequence", transform=ax.transAxes, color="#e07a7a", fontsize=12, va="bottom")
    _axfmt(ax)
    _save(fig, "sfms_all")


def fig_by_env(df, m, b):
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.4), sharex=True, sharey=True,
                             layout="constrained")
    # shared count scale across panels
    for k, ax in enumerate(axes.ravel()):
        sub = df[df["hard_class"] == k]
        draw_regions(ax, m, b)
        hb = ax.hexbin(sub["LOGMSTAR"], sub["logSFR"], gridsize=45, extent=[*MX, *MY],
                       mincnt=1, cmap="viridis", bins="log")
        ax.set_title(f"{CLASS_NAMES[k]}  (N={len(sub):,})",
                     color=COSMIC_WEB_COLORS[_CLASS_KEYS[k]], fontsize=15)
        ax.set_xlim(MX); ax.set_ylim(MY)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")
    _save(fig, "sfms_by_environment")


def fig_eigen(df, m, b):
    cols = ["lambda_1_mean", "lambda_2_mean", "lambda_3_mean"]
    labs = [r"$\langle\lambda_1\rangle$", r"$\langle\lambda_2\rangle$", r"$\langle\lambda_3\rangle$"]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), sharex=True, sharey=True,
                             layout="constrained")
    for ax, col, lab in zip(axes, cols, labs):
        draw_regions(ax, m, b, shade=False)
        # robust colour limits so sparse edge cells don't saturate the map
        vmin, vmax = np.percentile(df[col], [4, 96])
        hb = ax.hexbin(df["LOGMSTAR"], df["logSFR"], C=df[col], reduce_C_function=np.mean,
                       gridsize=38, extent=[*MX, *MY], mincnt=25, cmap="magma",
                       vmin=vmin, vmax=vmax)
        cb = fig.colorbar(hb, ax=ax); cb.set_label(f"mean {lab} per cell")
        ax.set_title(f"tidal environment: {lab}", fontsize=15)
        ax.set_xlim(MX); ax.set_ylim(MY)
        ax.set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
    axes[0].set_ylabel(r"$\log_{10}(\mathrm{SFR}/M_\odot\,\mathrm{yr}^{-1})$")
    _save(fig, "sfms_eigen_continuous")


def _save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=220, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", OUT / f"{stem}.*", flush=True)


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    df = load()
    m, b = fit_ms(df)
    print(f"MS fit: logSFR = {m:.3f} logM* + {b:.3f}   (N={len(df):,})", flush=True)
    fig_all(df, m, b)
    fig_by_env(df, m, b)
    fig_eigen(df, m, b)


if __name__ == "__main__":
    main()
