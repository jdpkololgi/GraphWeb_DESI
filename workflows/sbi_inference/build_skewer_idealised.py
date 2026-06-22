#!/usr/bin/env python3
"""Idealised fake-data skewer (slide-3 bookend) — the discrete→continuous intuition pump.

Same panel layout + standardised style as `build_skewer_animation.py` (it IMPORTS that
module's HTML template, so the bookend renders identically), but driven by a CLEAN
synthetic structure instead of real posteriors: the sightline crosses a single density
peak, so the three eigenvalue posteriors lift through λ_th in sequence
(void→wall→filament→cluster→…→void) and the class-probability bar morphs smoothly. No
real data — this is the conceptual setup; slide 10 then shows the *real* DESI version.

Eigenvalue means rise linearly with a density contrast d(s)∈[0,1] (Gaussian bump),
chosen so λ3, λ2, λ1 cross λ_th=0.2 at d≈0.33, 0.60, 0.86 — i.e. each extra eigenvalue
crossing the threshold flips the inferred class up one level (the ascending λ1≤λ2≤λ3
T-web convention used everywhere in the codebase).

Writes skewer_idealised.html into the SI run dir by default.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))
from build_skewer_animation import _HTML, CC, EC, ORDER, TH  # noqa: E402  reuse template + palette


def lam_means(d: float):
    """Ordered λ1≤λ2≤λ3 rising with density contrast d∈[0,1].

    Crossings of λ_th=0.2:  λ3 at d≈0.33 (→wall), λ2 at d≈0.60 (→filament),
    λ1 at d≈0.86 (→cluster). Below all crossings = void.
    """
    l3 = -0.10 + 0.90 * d
    l2 = -0.25 + 0.75 * d
    l1 = -0.40 + 0.70 * d
    return l1, l2, l3


def main(a):
    rng = np.random.default_rng(a.seed)
    smin, smax = 0.0, 200.0
    s0 = 0.5 * (smin + smax)
    centers = np.linspace(smin, smax, a.n_frames)
    xg = np.linspace(a.lmin, a.lmax, a.grid)

    def dens(s):  # density contrast along the skewer, peak 1.0 at the structure centre
        return np.exp(-0.5 * ((np.asarray(s) - s0) / a.struct_sigma) ** 2)

    frames = []
    for c in centers:
        d = float(dens(c))
        l1, l2, l3 = lam_means(d)
        sigp = a.post_sigma_min + (a.post_sigma_max - a.post_sigma_min) * d  # width grows where λ extreme
        kdes = [norm.pdf(xg, m, sigp).tolist() for m in (l1, l2, l3)]
        # Class probs from the ascending ordered posteriors (analytic, single Gaussians):
        p1 = float(1 - norm.cdf(TH, l1, sigp))
        p2 = float(1 - norm.cdf(TH, l2, sigp))
        p3 = float(1 - norm.cdf(TH, l3, sigp))
        probs = {"void": 1 - p3, "wall": p3 - p2, "filament": p2 - p1, "cluster": p1}
        frames.append({"s": float(c), "n": 50, "kde": kdes, "probs": probs, "width": float(sigp)})

    # Synthetic galaxies: denser near the peak; transverse scatter tightens in the clump;
    # colour = dominant inferred class (count of λ-means above λ_th).
    grid_s = np.linspace(smin, smax, 2000)
    w = 0.15 + dens(grid_s); w = w / w.sum()
    gs = rng.choice(grid_s, a.scatter_points, p=w)
    dd = dens(gs)
    tt = rng.normal(0.0, 6.0 * (1.0 - 0.6 * dd))
    cols = []
    for di in dd:
        l1, l2, l3 = lam_means(float(di))
        n_above = int(l1 > TH) + int(l2 > TH) + int(l3 > TH)
        cols.append(ORDER[n_above])
    gal = {"s": gs.round(2).tolist(), "t": tt.round(2).tolist(), "c": cols}

    payload = {"frames": frames, "xg": xg.tolist(), "gal": gal, "th": TH, "cc": CC, "ec": EC,
               "order": ORDER, "smin": float(centers[0]), "smax": float(centers[-1]),
               "tmin": float(np.percentile(tt, 2)), "tmax": float(np.percentile(tt, 98))}
    html = _HTML.replace("__DATA__", json.dumps(payload))
    html = html.replace("Line-of-sight skewer through the DESI wedge — real NPE posteriors",
                        "Idealised line-of-sight skewer — discrete → continuous concept")
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(html, encoding="utf-8")
    cl = [f["probs"]["cluster"] for f in frames]
    print(f"Saved: {out}  (peak P(cluster)={max(cl):.2f}; {a.n_frames} frames)", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
                    "desi_wedge_flowjax_linear_si/skewer_idealised.html")
    ap.add_argument("--n-frames", type=int, default=60)
    ap.add_argument("--struct-sigma", type=float, default=32.0, help="density-peak width (Mpc) along the skewer")
    ap.add_argument("--post-sigma-min", type=float, default=0.06, help="idealised posterior width in voids")
    ap.add_argument("--post-sigma-max", type=float, default=0.13, help="idealised posterior width in the cluster")
    ap.add_argument("--grid", type=int, default=120)
    ap.add_argument("--lmin", type=float, default=-1.0)
    ap.add_argument("--lmax", type=float, default=2.0)
    ap.add_argument("--scatter-points", type=int, default=6000)
    ap.add_argument("--seed", type=int, default=0)
    main(ap.parse_args())
