#!/usr/bin/env python3
"""Property->environment closure-test figures for the DESI LOA wedge (slide 10).

Reads the joined table from build_desi_wedge_property_join.py and makes themed
figures showing that known galaxy-evolution trends fall out of the NPE-inferred
cosmic-web environment -- truth-free evidence the posteriors carry real physical
information:

  Fig A  categorical 3-panel  : quenched fraction / median g-r / median log sSFR
                                vs the 4 inferred classes (void->cluster).
  Fig B  continuous  3-panel  : the same three properties vs a continuous
                                environment scalar E[Sum lambda] = trace of the
                                tidal tensor (~ local density); Spearman rho noted.
                                (--xaxis p_cluster for the P(cluster) version.)
  Fig C  mass control         : quenched fraction vs environment split into
                                stellar-mass bins (categorical + continuous), so
                                the trend cannot be dismissed as the mass-environment
                                relation.

Expected physics: quenched fraction and g-r RISE, log sSFR FALLS, toward
filaments/clusters and with increasing trace.

Theme: shared.plot_style (true-black, IBM Plex Sans, COSMIC_WEB_COLORS).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from shared.plot_style import (  # noqa: E402
    apply_style, finalize_axes, COSMIC_WEB_COLORS, CLASS_ORDER, ACCENT_COLORS,
)

ACCENT3 = [ACCENT_COLORS["blue"], ACCENT_COLORS["magenta"], ACCENT_COLORS["red"]]
_RNG = np.random.default_rng(0)


# --------------------------------------------------------------------------
# small stats helpers
# --------------------------------------------------------------------------
def wilson(k: int, n: int, z: float = 1.0):
    """Wilson score interval for a fraction (robust for small n / rare clusters)."""
    if n == 0:
        return np.nan, np.nan, np.nan
    p = k / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return p, centre - half, centre + half


def median_ci(v: np.ndarray, n_boot: int = 500, cap: int = 20000):
    """Bootstrap 1-sigma (16/84) CI on the median.

    The point estimate is taken on the same (optionally capped) subsample used
    for the bootstrap, and the interval is guarded to bracket it, so downstream
    asymmetric error whiskers are never negative.
    """
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan, np.nan
    s = v if v.size <= cap else _RNG.choice(v, cap, replace=False)
    med = float(np.median(s))
    if s.size < 3:
        return med, med, med
    boots = np.median(_RNG.choice(s, size=(n_boot, s.size), replace=True), axis=1)
    lo, hi = np.percentile(boots, [16, 84])
    return med, float(min(lo, med)), float(max(hi, med))


def yerr_pair(point: float, lo: float, hi: float):
    """Asymmetric, non-negative yerr for ax.errorbar (clamps tiny CI/point mismatches)."""
    return [[max(0.0, point - lo)], [max(0.0, hi - point)]]


def quantile_edges(x: np.ndarray, nbins: int = 10):
    x = x[np.isfinite(x)]
    edges = np.unique(np.quantile(x, np.linspace(0, 1, nbins + 1)))
    return edges


def binned_fraction(x: np.ndarray, flag: np.ndarray, edges: np.ndarray):
    idx = np.digitize(x, edges[1:-1])
    xc, p, lo, hi = [], [], [], []
    for b in range(len(edges) - 1):
        m = idx == b
        if m.sum() < 10:
            continue
        pp, l, h = wilson(int(flag[m].sum()), int(m.sum()))
        xc.append(np.median(x[m])); p.append(pp); lo.append(l); hi.append(h)
    return map(np.asarray, (xc, p, lo, hi))


def binned_median(x: np.ndarray, vals: np.ndarray, edges: np.ndarray):
    idx = np.digitize(x, edges[1:-1])
    xc, m_, lo, hi = [], [], [], []
    for b in range(len(edges) - 1):
        m = idx == b
        finite = m & np.isfinite(vals)
        if finite.sum() < 10:
            continue
        med, l, h = median_ci(vals[finite])
        xc.append(np.median(x[m])); m_.append(med); lo.append(l); hi.append(h)
    return map(np.asarray, (xc, m_, lo, hi))


def _rho(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 10:
        return np.nan
    return spearmanr(x[m], y[m]).correlation


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------
def fig_categorical(d, outdir: Path):
    """Fig A: property vs the 4 inferred classes."""
    hc = d["hard_class"]; xpos = np.arange(4)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    colors = [COSMIC_WEB_COLORS[c] for c in CLASS_ORDER]

    # quenched fraction
    for k in range(4):
        m = (hc == k) & d["has_props"]
        p, lo, hi = wilson(int(d["quenched"][m].sum()), int(m.sum()))
        axes[0].errorbar(k, p, yerr=yerr_pair(p, lo, hi), fmt="o", ms=10,
                         color=colors[k], capsize=4)
    finalize_axes(axes[0], "Quenched fraction vs inferred environment",
                  "inferred class", r"$f_{\rm quenched}$ (log sSFR $<-11$)", legend=False)

    # median g-r
    for k in range(4):
        m = (hc == k) & d["has_props"]
        med, lo, hi = median_ci(d["gr"][m])
        axes[1].errorbar(k, med, yerr=yerr_pair(med, lo, hi), fmt="s", ms=10,
                         color=colors[k], capsize=4)
    finalize_axes(axes[1], "Rest-frame colour vs inferred environment",
                  "inferred class", r"median $(g-r)_{0.1}$", legend=False)

    # median log sSFR
    for k in range(4):
        m = (hc == k) & d["has_ssfr"]
        med, lo, hi = median_ci(d["log_sSFR"][m])
        axes[2].errorbar(k, med, yerr=yerr_pair(med, lo, hi), fmt="D", ms=9,
                         color=colors[k], capsize=4)
    finalize_axes(axes[2], "Specific SFR vs inferred environment",
                  "inferred class", r"median $\log_{10}({\rm sSFR}\,/\,{\rm yr}^{-1})$", legend=False)

    for ax in axes:
        ax.set_xticks(xpos)
        ax.set_xticklabels([c.capitalize() for c in CLASS_ORDER])
    fig.suptitle("Closure test: galaxy properties vs inferred cosmic-web class (DESI LOA wedge)",
                 fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = outdir / "closure_categorical.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    print(f"[plot] wrote {p}", flush=True)


def fig_continuous(d, outdir: Path, xaxis: str):
    """Fig B: property vs a continuous environment scalar."""
    xname = {"trace": "trace_lambda", "p_cluster": "P_cluster"}[xaxis]
    xlabel = {"trace": r"$E[\,\lambda_1+\lambda_2+\lambda_3\,]$  (tidal-tensor trace $\propto$ density)",
              "p_cluster": r"$P(\mathrm{cluster})$"}[xaxis]
    x_all = d[xname]
    col = ACCENT_COLORS["magenta"]
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    mp = d["has_props"]
    edges = quantile_edges(x_all[mp], nbins=10)

    # quenched fraction
    xc, p, lo, hi = binned_fraction(x_all[mp], d["quenched"][mp], edges)
    axes[0].fill_between(xc, lo, hi, color=col, alpha=0.25)
    axes[0].plot(xc, p, "-o", color=col)
    rho = _rho(x_all[mp], d["quenched"][mp].astype(float))
    axes[0].text(0.05, 0.92, rf"$\rho_s={rho:+.2f}$", transform=axes[0].transAxes)
    finalize_axes(axes[0], "Quenched fraction vs environment", xlabel,
                  r"$f_{\rm quenched}$", legend=False)

    # median g-r
    xc, m_, lo, hi = binned_median(x_all[mp], d["gr"][mp], edges)
    axes[1].fill_between(xc, lo, hi, color=col, alpha=0.25)
    axes[1].plot(xc, m_, "-s", color=col)
    rho = _rho(x_all[mp], d["gr"][mp])
    axes[1].text(0.05, 0.92, rf"$\rho_s={rho:+.2f}$", transform=axes[1].transAxes)
    finalize_axes(axes[1], "Rest-frame colour vs environment", xlabel,
                  r"median $(g-r)_{0.1}$", legend=False)

    # median log sSFR
    ms = d["has_ssfr"]
    edges_s = quantile_edges(x_all[ms], nbins=10)
    xc, m_, lo, hi = binned_median(x_all[ms], d["log_sSFR"][ms], edges_s)
    axes[2].fill_between(xc, lo, hi, color=col, alpha=0.25)
    axes[2].plot(xc, m_, "-D", color=col)
    rho = _rho(x_all[ms], d["log_sSFR"][ms])
    axes[2].text(0.05, 0.08, rf"$\rho_s={rho:+.2f}$", transform=axes[2].transAxes)
    finalize_axes(axes[2], "Specific SFR vs environment", xlabel,
                  r"median $\log_{10}({\rm sSFR}\,/\,{\rm yr}^{-1})$", legend=False)

    fig.suptitle("Closure test: galaxy properties vs continuous inferred environment (DESI LOA wedge)",
                 fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = outdir / f"closure_continuous_{xaxis}.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    print(f"[plot] wrote {p}", flush=True)


def fig_mass_control(d, outdir: Path):
    """Fig C: quenched fraction vs environment in stellar-mass bins."""
    mp = d["has_props"]
    logm = d["LOGMSTAR"][mp]
    m_edges = np.quantile(logm, [0, 1 / 3, 2 / 3, 1.0])
    labels = [rf"$\log M_*\!\in[{m_edges[i]:.1f},{m_edges[i+1]:.1f})$" for i in range(3)]

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5))
    hc = d["hard_class"]; xpos = np.arange(4)

    # left: categorical, one line per mass bin
    for j in range(3):
        mb = mp & (d["LOGMSTAR"] >= m_edges[j]) & (d["LOGMSTAR"] < (m_edges[j + 1] + (1e-6 if j == 2 else 0)))
        ys, los, his = [], [], []
        for k in range(4):
            m = mb & (hc == k)
            p, lo, hi = wilson(int(d["quenched"][m].sum()), int(m.sum()))
            ys.append(p); los.append(lo); his.append(hi)
        ys, los, his = np.array(ys), np.array(los), np.array(his)
        axes[0].errorbar(xpos, ys, yerr=[np.maximum(0, ys - los), np.maximum(0, his - ys)],
                         fmt="-o", color=ACCENT3[j], capsize=3, label=labels[j])
    axes[0].set_xticks(xpos); axes[0].set_xticklabels([c.capitalize() for c in CLASS_ORDER])
    finalize_axes(axes[0], "Quenched fraction vs class, at fixed mass",
                  "inferred class", r"$f_{\rm quenched}$", legend=True)

    # right: continuous, one line per mass bin
    x_all = d["trace_lambda"]
    for j in range(3):
        mb = mp & (d["LOGMSTAR"] >= m_edges[j]) & (d["LOGMSTAR"] < (m_edges[j + 1] + (1e-6 if j == 2 else 0)))
        edges = quantile_edges(x_all[mb], nbins=6)
        xc, p, lo, hi = binned_fraction(x_all[mb], d["quenched"][mb], edges)
        axes[1].fill_between(xc, lo, hi, color=ACCENT3[j], alpha=0.18)
        axes[1].plot(xc, p, "-o", color=ACCENT3[j], label=labels[j])
    finalize_axes(axes[1], "Quenched fraction vs environment, at fixed mass",
                  r"$E[\,\lambda_1+\lambda_2+\lambda_3\,]$ ($\propto$ density)",
                  r"$f_{\rm quenched}$", legend=True)

    fig.suptitle("Mass control: the environment trend survives at fixed stellar mass",
                 fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = outdir / "closure_mass_control.png"
    fig.savefig(p, bbox_inches="tight"); plt.close(fig)
    print(f"[plot] wrote {p}", flush=True)


def load_table(path: Path):
    import pandas as pd
    if path.suffix == ".parquet" and not path.exists():
        path = path.with_suffix(".csv")
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
    # numpy views for the helpers
    cols = ["hard_class", "trace_lambda", "P_cluster", "gr", "log_sSFR",
            "LOGMSTAR", "quenched", "has_props", "has_ssfr"]
    return {c: df[c].to_numpy() for c in cols}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, default=None,
                    help="Joined parquet from the build step "
                         "(default: <preds dir>/desi_wedge_env_props.parquet next to the SI run).")
    ap.add_argument("--xaxis", choices=["trace", "p_cluster"], default="trace",
                    help="Continuous environment axis for Fig B (default: trace).")
    ap.add_argument("--outdir", type=Path, default=None,
                    help="Figure dir (default: <table dir>/closure/).")
    args = ap.parse_args()

    sys.path.insert(0, str(_REPO))
    from shared.config_paths import GRAPHWEB_SCRATCH_ROOT
    default_tbl = Path(f"{GRAPHWEB_SCRATCH_ROOT}/flowjax_inference_outputs/"
                       "desi_wedge_flowjax_linear_si/desi_wedge_env_props.parquet")
    table = (args.table or default_tbl).expanduser().resolve()
    outdir = (args.outdir or (table.parent / "closure")).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    apply_style()
    d = load_table(table)
    n_props = int(d["has_props"].sum())
    print(f"[plot] table {table}  ({d['hard_class'].size:,} rows, {n_props:,} with valid mass)", flush=True)

    fig_categorical(d, outdir)
    fig_continuous(d, outdir, args.xaxis)
    fig_mass_control(d, outdir)
    print(f"[plot] done -> {outdir}", flush=True)


if __name__ == "__main__":
    main()
