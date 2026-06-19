#!/usr/bin/env python3
"""Themed figures + truth-free validations for the DESI FlowJAX NPE inference.

Reads desi_wedge_flowjax_preds.npz + summary.json from infer_desi_wedge_flowjax.py
(no GPU). Produces, in the plot-style-guide theme (true-black, COSMIC_WEB_COLORS):
  1. class_fractions_comparison.png  — DESI NPE vs Abacus truth vs regression DESI
  2. class_sky_map.png               — RA/Dec by hard class + soft P(filament)
  3. posterior_width_sky_map.png     — RA/Dec by posterior width + entropy (NPE-only)
  4. domain_shift_overlay.png        — DESI inferred λ vs Abacus training λ
  5. width_vs_boundary.png           — posterior width vs distance to survey edge
                                       and to the T-web class boundary (λ≈λ_th)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILLUSTRIS_ROOT = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILLUSTRIS_ROOT) not in sys.path:
    sys.path.insert(0, str(_ILLUSTRIS_ROOT))
from shared.plot_style import apply_style, COSMIC_WEB_COLORS, CLASS_ORDER, ACCENT_COLORS  # noqa: E402

CLASS_CMAP = matplotlib.colors.ListedColormap([COSMIC_WEB_COLORS[c] for c in CLASS_ORDER])


def _entropy(p):
    p = np.clip(p, 1e-12, 1.0)
    return -(p * np.log(p)).sum(axis=1)


def main(args):
    apply_style()
    d = np.load(args.preds_npz)
    summary = json.loads(Path(args.summary_json).read_text())
    out = Path(args.output_dir or Path(args.preds_npz).parent)
    out.mkdir(parents=True, exist_ok=True)
    lam_th = float(summary.get("lambda_threshold", 0.2))

    ra, dec, z = d["ra"], d["dec"], d["z"]
    lam_mean, lam_std = d["lambda_mean"], d["lambda_std"]   # [N,3]
    classprob, hard = d["classprob"], d["hard_class"]       # [N,4],[N]
    width = lam_std.mean(axis=1)                            # scalar posterior width
    ent = _entropy(classprob)

    # ---- 1. class-fraction comparison (4 series: sim truth/sim NPE/real reg/real NPE) ----
    npe = summary["desi"]["class_fractions_npe"]
    ref = summary.get("reference_fractions", {})
    truth = ref.get("abacus_truth") or {}
    reg = ref.get("regression_desi") or {}
    # Abacus NPE (the trained model's predictions ON the Abacus test set) — hard argmax
    # fractions from the saved per-galaxy class-prob npz. This 4th series localizes the
    # cluster gap: sim-truth ≈ sim-NPE but real-reg ≈ real-NPE => gap enters at transfer.
    abacus_npe = {}
    if args.abacus_classprob_npz:
        an = np.load(args.abacus_classprob_npz)
        P = np.stack([np.asarray(an[c]) for c in CLASS_ORDER], axis=1)
        h = P.argmax(axis=1)
        abacus_npe = {c: float(np.mean(h == k)) for k, c in enumerate(CLASS_ORDER)}
    series = [("Abacus truth", truth, 0.4), ("Abacus NPE", abacus_npe, 0.62),
              ("Regression DESI", reg, 0.82), ("NPE DESI", npe, 1.0)]
    fig, ax = plt.subplots(figsize=(11, 6))
    x = np.arange(4); w = 0.2
    cols = [COSMIC_WEB_COLORS[c] for c in CLASS_ORDER]
    from matplotlib.patches import Patch
    for i, (lab, frac, alpha) in enumerate(series):
        vals = [float(frac.get(c, np.nan)) for c in CLASS_ORDER]
        bars = ax.bar(x + (i - 1.5) * w, vals, w, color=cols, alpha=alpha, edgecolor="#F2F2F2", lw=0.7)
        for b, v in zip(bars, vals):
            if np.isfinite(v):
                ax.text(b.get_x() + b.get_width() / 2, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels([c.capitalize() for c in CLASS_ORDER])
    ax.set_ylabel("Class fraction"); ax.set_ylim(0.0, 1.0)
    ax.set_title(rf"T-Web class fractions: simulation vs real DESI ($\lambda_{{th}}$={lam_th})")
    # Bar COLOR = class (void/wall/filament/cluster); bar SHADE = method. The legend
    # must encode the shade, so use neutral swatches at the four alpha levels.
    handles = [Patch(facecolor="#F2F2F2", edgecolor="#F2F2F2", alpha=a, label=lab) for (lab, _, a) in series]
    ax.legend(handles=handles, ncol=2, title="shade = method  (colour = class)")
    fig.savefig(out / "class_fractions_comparison.png", bbox_inches="tight"); plt.close(fig)

    # ---- 2. class sky map + soft P(filament) ----
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    sc = axes[0].scatter(ra, dec, c=hard, cmap=CLASS_CMAP, s=2, vmin=-0.5, vmax=3.5, rasterized=True)
    axes[0].set_title("Inferred T-Web class"); axes[0].set_xlabel("RA [deg]"); axes[0].set_ylabel("Dec [deg]")
    cb = fig.colorbar(sc, ax=axes[0], ticks=[0, 1, 2, 3]); cb.ax.set_yticklabels([c.capitalize() for c in CLASS_ORDER])
    p_fil = classprob[:, CLASS_ORDER.index("filament")]
    sc2 = axes[1].scatter(ra, dec, c=p_fil, cmap="magma", s=2, vmin=0, vmax=1, rasterized=True)
    axes[1].set_title("Posterior P(filament)"); axes[1].set_xlabel("RA [deg]"); axes[1].set_ylabel("Dec [deg]")
    fig.colorbar(sc2, ax=axes[1], label="P(filament)")
    fig.savefig(out / "class_sky_map.png", bbox_inches="tight", dpi=200); plt.close(fig)

    # ---- 3. posterior-width / uncertainty sky map (the NPE advantage) ----
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    sc = axes[0].scatter(ra, dec, c=width, cmap="viridis", s=2,
                         vmin=np.percentile(width, 2), vmax=np.percentile(width, 98), rasterized=True)
    axes[0].set_title("Posterior width (mean λ std)"); axes[0].set_xlabel("RA [deg]"); axes[0].set_ylabel("Dec [deg]")
    fig.colorbar(sc, ax=axes[0], label=r"mean $\sigma_\lambda$")
    sc2 = axes[1].scatter(ra, dec, c=ent, cmap="viridis", s=2,
                          vmin=np.percentile(ent, 2), vmax=np.percentile(ent, 98), rasterized=True)
    axes[1].set_title("Class-probability entropy"); axes[1].set_xlabel("RA [deg]"); axes[1].set_ylabel("Dec [deg]")
    fig.colorbar(sc2, ax=axes[1], label="entropy [nats]")
    fig.savefig(out / "posterior_width_sky_map.png", bbox_inches="tight", dpi=200); plt.close(fig)

    # ---- 4. domain-shift overlay: DESI inferred λ vs Abacus training λ ----
    abacus_eig = None
    if args.calibration_cache:
        import pickle
        with open(args.calibration_cache, "rb") as f:
            abacus_eig = np.asarray(pickle.load(f).get("eigenvalues_raw"))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    accent = [ACCENT_COLORS["blue"], ACCENT_COLORS["magenta"], ACCENT_COLORS["red"]]
    for k in range(3):
        axes[k].hist(lam_mean[:, k], bins=80, density=True, alpha=0.8, color=accent[k], label="DESI NPE (post. mean)")
        if abacus_eig is not None:
            axes[k].hist(abacus_eig[:, k], bins=80, density=True, histtype="step", lw=2, color="#F2F2F2", label="Abacus train")
        axes[k].axvline(lam_th, ls=":", color="#F2F2F2", alpha=0.6)
        axes[k].set_xlabel(rf"$\lambda_{k+1}$"); axes[k].set_ylabel("density" if k == 0 else "")
        if k == 0:
            axes[k].legend(fontsize=9)
    fig.suptitle("Domain-shift check: DESI inferred λ vs Abacus training λ")
    fig.savefig(out / "domain_shift_overlay.png", bbox_inches="tight"); plt.close(fig)

    # ---- 5. width vs boundary distance (survey edge + class boundary) ----
    # normalized distance to nearest wedge footprint edge (RA/Dec/z), in [0, .5]
    def norm_edge_dist(v, lo, hi):
        return np.minimum(v - lo, hi - v) / (hi - lo)
    edge_d = np.minimum.reduce([
        norm_edge_dist(ra, args.ra_min, args.ra_max),
        norm_edge_dist(dec, args.dec_min, args.dec_max),
        norm_edge_dist(z, args.z_min, args.z_max)])
    class_d = np.abs(lam_mean - lam_th).min(axis=1)  # distance to T-web class boundary

    def binned(xv, yv, nb=25):
        edges = np.quantile(xv, np.linspace(0, 1, nb + 1))
        idx = np.clip(np.digitize(xv, edges[1:-1]), 0, nb - 1)
        cx = np.array([xv[idx == b].mean() if np.any(idx == b) else np.nan for b in range(nb)])
        cy = np.array([yv[idx == b].mean() if np.any(idx == b) else np.nan for b in range(nb)])
        return cx, cy

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bx, by = binned(edge_d, width)
    axes[0].plot(bx, by, "-o", color=ACCENT_COLORS["magenta"], ms=4)
    axes[0].set_xlabel("normalized distance to survey edge"); axes[0].set_ylabel(r"mean posterior width $\sigma_\lambda$")
    axes[0].set_title("Width vs survey-edge distance")
    bx2, by2 = binned(class_d, width)
    axes[1].plot(bx2, by2, "-o", color=ACCENT_COLORS["blue"], ms=4)
    axes[1].set_xlabel(r"distance to class boundary $\min_k|\lambda_k-\lambda_{th}|$")
    axes[1].set_ylabel(r"mean posterior width $\sigma_\lambda$")
    axes[1].set_title("Width vs class-boundary distance")
    fig.savefig(out / "width_vs_boundary.png", bbox_inches="tight"); plt.close(fig)

    # ---- 6. 3D interactive maps (plotly) — clearer than the dense 2D sky scatter ----
    try:
        import plotly.graph_objects as go
        from astropy.cosmology import Planck18 as cosmo
        dist = cosmo.comoving_distance(z).value
        rar, decr = np.deg2rad(ra), np.deg2rad(dec)
        X = dist * np.cos(decr) * np.cos(rar); Y = dist * np.cos(decr) * np.sin(rar); Zc = dist * np.sin(decr)
        # downsample for a responsive HTML
        n = len(X); sel = np.random.default_rng(0).choice(n, min(n, args.max_3d_points), replace=False)
        # (a) by class
        fig3 = go.Figure()
        for k, c in enumerate(CLASS_ORDER):
            m = sel[hard[sel] == k]
            fig3.add_trace(go.Scatter3d(x=X[m], y=Y[m], z=Zc[m], mode="markers",
                           marker=dict(size=1.6, color=COSMIC_WEB_COLORS[c]), name=c.capitalize()))
        fig3.update_layout(template="plotly_dark", title="DESI wedge — inferred T-Web class (NPE)",
                           scene=dict(xaxis_title="X [Mpc]", yaxis_title="Y [Mpc]", zaxis_title="Z [Mpc]"),
                           paper_bgcolor="#000000")
        fig3.write_html(str(out / "class_3d.html"), include_plotlyjs="cdn")
        # (b) by posterior width
        fig4 = go.Figure(go.Scatter3d(
            x=X[sel], y=Y[sel], z=Zc[sel], mode="markers",
            marker=dict(size=1.6, color=width[sel], colorscale="Viridis", colorbar=dict(title="σ_λ"),
                        cmin=np.percentile(width, 2), cmax=np.percentile(width, 98))))
        fig4.update_layout(template="plotly_dark", title="DESI wedge — posterior width (NPE uncertainty)",
                           scene=dict(xaxis_title="X [Mpc]", yaxis_title="Y [Mpc]", zaxis_title="Z [Mpc]"),
                           paper_bgcolor="#000000")
        fig4.write_html(str(out / "posterior_width_3d.html"), include_plotlyjs="cdn")
        print("  3D:", out / "class_3d.html", "+ posterior_width_3d.html")
    except ImportError as e:
        print(f"  (plotly/astropy unavailable, skipping 3D HTML: {e})")

    print("Saved figures to", out)
    for p in ["class_fractions_comparison", "class_sky_map", "posterior_width_sky_map",
              "domain_shift_overlay", "width_vs_boundary"]:
        print("  ", out / f"{p}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Themed DESI FlowJAX NPE figures")
    ap.add_argument("--preds-npz", required=True)
    ap.add_argument("--summary-json", required=True)
    ap.add_argument("--calibration-cache", default=None, help="for the Abacus-λ domain-shift overlay")
    ap.add_argument("--abacus-classprob-npz", default=None,
                    help="Abacus NPE per-galaxy class probs (4th class-fraction series)")
    ap.add_argument("--output-dir", default=None)
    ap.add_argument("--ra-min", type=float, default=120.0); ap.add_argument("--ra-max", type=float, default=160.0)
    ap.add_argument("--dec-min", type=float, default=14.5); ap.add_argument("--dec-max", type=float, default=30.6)
    ap.add_argument("--z-min", type=float, default=0.2); ap.add_argument("--z-max", type=float, default=0.3)
    ap.add_argument("--max-3d-points", type=int, default=60000)
    main(ap.parse_args())
