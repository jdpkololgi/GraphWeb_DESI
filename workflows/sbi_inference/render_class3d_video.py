#!/usr/bin/env python3
"""Rotating fly-through GIF of the inferred DESI cosmic-web environments (3D).

The interactive class_3d.html is a plotly canvas Keynote can't play; this rebuilds the
same point cloud — comoving XYZ from preds.npz (ra/dec/z), coloured by inferred T-web
class — and renders a matplotlib 3D scatter with a slow rotating + gently zooming camera,
saved as a looping GIF (Keynote-native).

Build-up: the environments fade in successively (void -> wall -> filament -> cluster), so
the cosmic web assembles from diffuse to dense before the full rotation. All clusters are
kept (rare class); the other classes are subsampled for speed.

Writes class_3d_flythrough.gif into the SI run dir by default.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402
from astropy.cosmology import Planck18 as cosmo  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from shared.plot_style import apply_style, COSMIC_WEB_COLORS, CLASS_ORDER  # noqa: E402

DEF_PREDS = ("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
             "desi_wedge_flowjax_linear_si/desi_wedge_flowjax_preds.npz")

# proportionate sizes / opacities (cluster only slightly larger so the rare knots read)
SIZES = {0: 1.8, 1: 2.0, 2: 2.4, 3: 3.2}
ALPHAS = {0: 0.30, 1: 0.42, 2: 0.55, 3: 0.90}
ADD_LABEL = ["Voids", "+ Walls", "+ Filaments", "+ Clusters"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--preds", type=Path, default=Path(DEF_PREDS))
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--n-points", type=int, default=28000, help="subsample of the non-cluster classes")
    ap.add_argument("--n-frames", type=int, default=96)
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--reveal-frac", type=float, default=0.6, help="fraction of frames for the build-up (0=off)")
    ap.add_argument("--rotations", type=float, default=1.2, help="total camera revolutions")
    ap.add_argument("--dpi", type=int, default=110)
    ap.add_argument("--fig-w", type=float, default=12.0, help="figure width inches (widescreen)")
    ap.add_argument("--fig-h", type=float, default=6.75, help="figure height inches (16:9)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()
    out = args.out or (args.preds.parent / "class_3d_flythrough.gif")

    d = np.load(args.preds)
    ra, dec, z = (np.asarray(d[k], float) for k in ("ra", "dec", "z"))
    hard = np.asarray(d["hard_class"])
    dist = cosmo.comoving_distance(z).value
    rar, decr = np.radians(ra), np.radians(dec)
    X = dist * np.cos(decr) * np.cos(rar)
    Y = dist * np.cos(decr) * np.sin(rar)
    Z = dist * np.sin(decr)

    rng = np.random.default_rng(args.seed)               # keep all clusters; subsample the rest
    rest = np.where(hard != 3)[0]
    sel = np.concatenate([rng.choice(rest, min(len(rest), args.n_points), replace=False),
                          np.where(hard == 3)[0]])
    Xs, Ys, Zs, hs = X[sel], Y[sel], Z[sel], hard[sel]

    apply_style()
    fig = plt.figure(figsize=(args.fig_w, args.fig_h)); fig.patch.set_facecolor("#000000")
    ax = fig.add_subplot(111, projection="3d"); ax.set_facecolor("#000000")
    colls = []
    for k, c in enumerate(CLASS_ORDER):           # void..cluster so clusters draw on top
        m = hs == k
        colls.append(ax.scatter(Xs[m], Ys[m], Zs[m], s=SIZES[k], c=COSMIC_WEB_COLORS[c],
                                alpha=ALPHAS[k], edgecolors="none", depthshade=False))

    cx, cy, cz = Xs.mean(), Ys.mean(), Zs.mean()
    R = 0.5 * max(np.ptp(Xs), np.ptp(Ys), np.ptp(Zs))
    ax.set_box_aspect((1, 1, 1), zoom=1.4)    # enlarge cube (slab corners are empty) to fill the frame
    ax.set_position([0.0, 0.06, 1.0, 0.83])   # full width; reserve top (title+caption) and bottom (legend) bands
    ax.set_axis_off()
    handles = [Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=COSMIC_WEB_COLORS[c],
                      markeredgecolor="none", markersize=10, label=c.capitalize()) for c in CLASS_ORDER]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.012), ncol=4, fontsize=12,
               frameon=False, labelcolor="#F2F2F2", columnspacing=1.8, handletextpad=0.4)
    fig.suptitle("DESI BGS wedge — inferred cosmic-web environment", color="#F2F2F2", fontsize=15, y=0.965)
    cap = fig.text(0.5, 0.905, "", ha="center", va="top", color="#F2F2F2", fontsize=13)

    n = args.n_frames
    reveal = int(n * args.reveal_frac)
    per = max(1, reveal // 4)
    fade = max(1.0, per * 0.6)

    def class_alpha(k, i):
        if reveal <= 0:
            return ALPHAS[k]
        start = k * per
        if i < start:
            return 0.0
        return ALPHAS[k] * min(1.0, (i - start) / fade)

    def caption(i):
        if reveal <= 0 or i >= reveal:
            return "Inferred cosmic web"
        return ADD_LABEL[min(3, i // per)]

    def update(i):
        t = i / n
        for k, coll in enumerate(colls):
            coll.set_alpha(class_alpha(k, i))
        cap.set_text(caption(i))
        ax.view_init(elev=18.0 + 10.0 * np.sin(2 * np.pi * t),
                     azim=360.0 * args.rotations * t)
        s = 0.88 + 0.12 * np.cos(2 * np.pi * t)          # gentle fly in/out
        ax.set_xlim(cx - R * s, cx + R * s); ax.set_ylim(cy - R * s, cy + R * s)
        ax.set_zlim(cz - R * s, cz + R * s)
        return ()

    ani = FuncAnimation(fig, update, frames=n, interval=1000 / args.fps, blit=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    ani.save(str(out), writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig)
    print(f"[class3d] wrote {out}  ({n} frames @ {args.fps} fps, {len(Xs):,} points, "
          f"{int((hs==3).sum()):,} clusters, reveal={args.reveal_frac})", flush=True)


if __name__ == "__main__":
    main()
