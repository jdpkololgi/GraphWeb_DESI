#!/usr/bin/env python3
"""Sweep the T-web threshold λ_th and show how the four class fractions migrate
(cluster->filament->wall->void as λ_th rises), comparing Abacus truth, Abacus NPE
and DESI NPE. Outputs a static line plot and an interactive animation.
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, COSMIC_WEB_COLORS, CLASS_ORDER  # noqa: E402


def fractions(X, th):
    n = (X > th).sum(1)            # ascending eigenvalues -> n = number above
    return np.array([np.mean(n == k) for k in range(4)])


def main(a):
    apply_style()
    aself = np.load(a.abacus_self_npz); desi = np.load(a.desi_npz)
    series = [("Abacus truth", aself["eig_truth"], "-"),
              ("Abacus NPE", aself["lambda_mean"], "--"),
              ("DESI NPE", desi["lambda_mean"], ":")]
    ths = np.round(np.linspace(a.th_min, a.th_max, a.n_th), 3)
    outdir = Path(a.desi_npz).parent

    # precompute fractions[series][th][class]
    data = {lab: np.array([fractions(X, t) for t in ths]) for lab, X, _ in series}

    # ---- static line plot: 1 panel per class, 3 lines ----
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.5), sharey=True)
    for k, c in enumerate(CLASS_ORDER):
        for lab, _, ls in series:
            axes[k].plot(ths, data[lab][:, k], ls, lw=2, color=COSMIC_WEB_COLORS[c], label=lab)
        axes[k].axvline(0.2, color="#888", ls=":", lw=1)
        axes[k].set_title(c.capitalize()); axes[k].set_xlabel(r"$\lambda_{th}$")
        if k == 0:
            axes[k].set_ylabel("class fraction"); axes[k].legend(fontsize=8)
    fig.suptitle("T-web class fractions vs threshold  (line = class colour; style = method)")
    fig.savefig(outdir / "lambda_th_sweep.png", bbox_inches="tight", dpi=200); plt.close(fig)
    print("Saved:", outdir / "lambda_th_sweep.png")

    # ---- interactive animation: grouped bar morphing with λ_th ----
    try:
        import plotly.graph_objects as go
        labels = [s[0] for s in series]
        alphas = [0.45, 0.7, 1.0]
        def bars(ti):
            return [go.Bar(name=lab, x=[c.capitalize() for c in CLASS_ORDER],
                           y=list(data[lab][ti]),
                           marker=dict(color=[COSMIC_WEB_COLORS[c] for c in CLASS_ORDER], opacity=al))
                    for (lab, _, _), al in zip(series, alphas)]
        frames = [go.Frame(data=bars(ti), name=f"{ths[ti]:.2f}") for ti in range(len(ths))]
        i0 = int(np.argmin(np.abs(ths - 0.2)))
        fig = go.Figure(data=bars(i0), frames=frames)
        fig.update_layout(
            template="plotly_dark", paper_bgcolor="#000000", barmode="group",
            title="T-web class fractions as the threshold λ_th sweeps",
            yaxis=dict(title="class fraction", range=[0, 1]),
            updatemenus=[dict(type="buttons", showactive=False, x=0.0, y=1.15,
                              buttons=[dict(label="play", method="animate",
                                            args=[None, dict(frame=dict(duration=120, redraw=True), fromcurrent=True)]),
                                       dict(label="pause", method="animate",
                                            args=[[None], dict(frame=dict(duration=0), mode="immediate")])])],
            sliders=[dict(active=i0, currentvalue=dict(prefix="λ_th = "),
                          steps=[dict(method="animate", label=f"{ths[ti]:.2f}",
                                      args=[[f"{ths[ti]:.2f}"], dict(mode="immediate", frame=dict(duration=0, redraw=True))])
                                 for ti in range(len(ths))])])
        fig.write_html(str(outdir / "lambda_th_sweep_animation.html"), include_plotlyjs="cdn")
        print("Saved:", outdir / "lambda_th_sweep_animation.html")
    except ImportError as e:
        print("plotly unavailable:", e)

    # report cluster fraction at a couple of thresholds
    for t in (0.1, 0.2, 0.3):
        ti = int(np.argmin(np.abs(ths - t)))
        print(f"  λ_th={t}: cluster frac  " + "  ".join(f"{lab}={data[lab][ti,3]:.3f}" for lab in [s[0] for s in series]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--abacus-self-npz", required=True)
    ap.add_argument("--desi-npz", required=True)
    ap.add_argument("--th-min", type=float, default=-0.1); ap.add_argument("--th-max", type=float, default=0.6)
    ap.add_argument("--n-th", type=int, default=36)
    main(ap.parse_args())
