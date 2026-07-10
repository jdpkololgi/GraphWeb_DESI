#!/usr/bin/env python3
"""Continuous / fixed-mass views of environment vs galaxy properties (CIGALE-HZ).

1. env_mass_surface_3d      3D surface: quenched fraction over (logM*, tidal trace)
2. env_mass_heatmaps        2D heatmaps: median (g-r)/log sSFR/quenched over (logM*, trace)
3. env_class_by_massbin     static: property vs 4 classes, one line per stellar-mass bin
4. env_class_mass_animation animated mass slider: 4 environments' g-r/sSFR/SFR vs rising M*
"""
import os, sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.lines import Line2D

sys.path.insert(0, "/global/homes/d/dkololgi/TNG/Illustris")
from shared.plot_style import apply_style, COSMIC_WEB_COLORS, ACCENT_COLORS  # noqa

PARQUET = Path(os.environ.get("WEDGE_PARQUET",
    "/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
    "desi_wedge_cigale_hz/desi_wedge_env_props.parquet"))
OUT = Path(os.environ.get("WEDGE_FIGDIR",
    "/pscratch/sd/d/dkololgi/graphweb_desi/figures/desi_wedge_cigale_hz"))
CKEYS = ["void", "wall", "filament", "cluster"]; CNAMES = ["Void", "Wall", "Filament", "Cluster"]
CCOL = [COSMIC_WEB_COLORS[k] for k in CKEYS]
MLO, MHI = 10.4, 11.5     # CIGALE stellar-mass working range


def load():
    d = pd.read_parquet(PARQUET)
    d = d[d["has_ssfr"] & np.isfinite(d["LOGMSTAR"]) & np.isfinite(d["log_sSFR"])].copy()
    d["logSFR"] = np.log10(d["SFR"])
    return d


def _grid_stat(x, y, z, xe, ye, how, mincnt=30):
    ix = np.clip(np.digitize(x, xe) - 1, 0, len(xe) - 2)
    iy = np.clip(np.digitize(y, ye) - 1, 0, len(ye) - 2)
    G = np.full((len(ye) - 1, len(xe) - 1), np.nan)
    for i in range(len(ye) - 1):
        for j in range(len(xe) - 1):
            m = (ix == j) & (iy == i)
            if m.sum() >= mincnt:
                G[i, j] = z[m].mean() if how == "frac" else np.nanmedian(z[m])
    return G


def surface_3d(d):
    from mpl_toolkits.mplot3d import Axes3D  # noqa
    xe = np.linspace(MLO, MHI, 13)
    ye = np.quantile(d["trace_lambda"], np.linspace(0.02, 0.98, 13))
    G = _grid_stat(d["LOGMSTAR"].to_numpy(), d["trace_lambda"].to_numpy(),
                   d["quenched"].to_numpy().astype(float), xe, ye, "frac")
    xc = 0.5 * (xe[:-1] + xe[1:]); yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(xc, yc)
    # fill sparse holes by column-nan-interp for a smooth surface
    for i in range(G.shape[0]):
        row = G[i]; ok = np.isfinite(row)
        if ok.sum() >= 2:
            G[i] = np.interp(xc, xc[ok], row[ok])
    fig = plt.figure(figsize=(9.5, 7.5))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(X, Y, G, cmap="magma", vmin=np.nanmin(G), vmax=np.nanmax(G),
                           rstride=1, cstride=1, antialiased=True, alpha=0.96)
    ax.set_xlabel("\n" + r"$\log_{10}(M_*/M_\odot)$", fontsize=12)
    ax.set_ylabel("\n" + r"tidal trace  $E[\Sigma\lambda]$", fontsize=12)
    ax.set_zlabel("\nquenched fraction", fontsize=12)
    ax.set_title("Quenching over mass × environment: strong along mass,\n"
                 "gentle-but-monotonic along tidal trace", fontsize=13)
    ax.view_init(elev=26, azim=-60)
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.set_pane_color((0, 0, 0, 0))
    fig.colorbar(surf, ax=ax, shrink=0.5, pad=0.08, label="quenched fraction")
    _save(fig, "env_mass_surface_3d")


def heatmaps(d):
    xe = np.linspace(MLO, MHI, 15)
    ye = np.quantile(d["trace_lambda"], np.linspace(0.02, 0.98, 15))
    specs = [("gr", "med", r"median $(g-r)$", "viridis"),
             ("log_sSFR", "med", r"median $\log_{10}$ sSFR", "viridis_r"),
             ("quenched", "frac", "quenched fraction", "magma")]
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8), layout="constrained")
    for ax, (col, how, lab, cm) in zip(axes, specs):
        z = d[col].to_numpy().astype(float)
        G = _grid_stat(d["LOGMSTAR"].to_numpy(), d["trace_lambda"].to_numpy(), z, xe, ye, how)
        im = ax.imshow(G, origin="lower", aspect="auto", cmap=cm,
                       extent=[xe[0], xe[-1], ye[0], ye[-1]])
        fig.colorbar(im, ax=ax, label=lab)
        ax.set_xlabel(r"$\log_{10}(M_*/M_\odot)$")
        ax.set_ylabel(r"tidal trace  $E[\Sigma\lambda]$")
        ax.set_title(lab, fontsize=13)
    fig.suptitle("Galaxy properties over stellar mass × inferred environment (continuous)", fontsize=14)
    _save(fig, "env_mass_heatmaps")


def _class_medians(d, prop, how, msel):
    out, err = [], []
    for k in range(4):
        v = d.loc[msel & (d["hard_class"] == k), prop].to_numpy()
        v = v[np.isfinite(v)]
        if len(v) < 15:
            out.append(np.nan); err.append(0); continue
        if how == "frac":
            p = v.mean(); out.append(p); err.append(np.sqrt(p*(1-p)/len(v)))
        else:
            out.append(np.median(v))
            err.append(1.253 * np.std(v) / np.sqrt(len(v)))
    return np.array(out), np.array(err)


def static_massbins(d):
    edges = np.array([10.4, 10.7, 10.9, 11.1, 11.5])
    cmap = plt.cm.plasma(np.linspace(0.15, 0.9, len(edges) - 1))
    specs = [("gr", "med", r"median $(g-r)$"),
             ("log_sSFR", "med", r"median $\log_{10}$ sSFR"),
             ("logSFR", "med", r"median $\log_{10}$ SFR")]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), layout="constrained")
    xx = np.arange(4)
    for ax, (col, how, lab) in zip(axes, specs):
        for b in range(len(edges) - 1):
            msel = (d["LOGMSTAR"] >= edges[b]) & (d["LOGMSTAR"] < edges[b+1])
            y, e = _class_medians(d, col, how, msel)
            ax.errorbar(xx, y, yerr=e, marker="o", ms=6, lw=2, color=cmap[b],
                        label=rf"${edges[b]:.1f}\!-\!{edges[b+1]:.1f}$")
        ax.set_xticks(xx); ax.set_xticklabels(CNAMES, fontsize=11)
        ax.set_ylabel(lab, fontsize=12); ax.set_title(f"{lab} vs environment, per $\\log M_*$ bin", fontsize=12)
    axes[0].legend(title=r"$\log_{10}M_*$", fontsize=9, title_fontsize=9, frameon=False)
    _save(fig, "env_class_by_massbin")


def animate_massbins(d):
    centers = np.linspace(10.55, 11.35, 22); half = 0.15
    specs = [("gr", "med", r"median $(g-r)$"),
             ("log_sSFR", "med", r"median $\log_{10}$ sSFR"),
             ("logSFR", "med", r"median $\log_{10}$ SFR")]
    # fixed y-limits from global percentiles for stable framing
    ylims = []
    for col, _, _ in specs:
        v = d[col].to_numpy(); v = v[np.isfinite(v)]
        ylims.append((np.percentile(v, 2), np.percentile(v, 98)))
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2))
    fig.subplots_adjust(top=0.82, wspace=0.32, left=0.06, right=0.98, bottom=0.12)
    sup = fig.suptitle("", fontsize=15)
    xx = np.arange(4)
    arts = []
    for ax, (col, how, lab), yl in zip(axes, specs, ylims):
        ax.set_xticks(xx); ax.set_xticklabels(CNAMES, fontsize=11)
        ax.set_ylim(*yl); ax.set_ylabel(lab, fontsize=12)
        line, = ax.plot(xx, np.full(4, np.nan), "-", lw=2, color="#888")
        pts = ax.scatter(xx, np.full(4, np.nan), s=120, c=CCOL, zorder=5, edgecolors="k", linewidths=0.6)
        arts.append((ax, col, how, line, pts))

    def update(i):
        c = centers[i]
        msel = (d["LOGMSTAR"] >= c - half) & (d["LOGMSTAR"] < c + half)
        sup.set_text(rf"$\log_{{10}}(M_*/M_\odot) = {c:.2f}$   "
                     rf"(N={int(msel.sum()):,})   —   environment at fixed mass")
        for ax, col, how, line, pts in arts:
            y, _ = _class_medians(d, col, how, msel)
            line.set_ydata(y); pts.set_offsets(np.c_[xx, y])
        return [a[3] for a in arts] + [a[4] for a in arts] + [sup]

    apply_style()  # re-apply after fig creation for dark bg
    fig.set_facecolor("black")
    for ax, *_ in [(a[0],) for a in arts]:
        pass
    anim = animation.FuncAnimation(fig, update, frames=len(centers), blit=False)
    gif = OUT / "env_class_mass_animation.gif"
    anim.save(str(gif), writer=animation.PillowWriter(fps=5), dpi=90,
              savefig_kwargs={"facecolor": "black"})
    plt.close(fig); print("Saved:", gif)


def _save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{stem}.{ext}", dpi=200, bbox_inches="tight", facecolor="black")
    plt.close(fig); print("Saved:", OUT / f"{stem}.png")


def main():
    apply_style(); OUT.mkdir(parents=True, exist_ok=True)
    d = load(); print(f"N={len(d):,}  logM* {d['LOGMSTAR'].min():.2f}-{d['LOGMSTAR'].max():.2f}")
    surface_3d(d); heatmaps(d); static_massbins(d); animate_massbins(d)


if __name__ == "__main__":
    main()
