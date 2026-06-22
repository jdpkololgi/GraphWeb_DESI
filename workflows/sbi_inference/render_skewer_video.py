#!/usr/bin/env python3
"""Render a skewer HTML animation to mp4/gif for Keynote (no browser needed).

Keynote can't play the interactive <canvas> HTML, but it embeds mp4 (plays inline)
or gif. Each skewer HTML already carries its full data payload as embedded JSON
(`var D={...}`), so we extract that and re-render the SAME panels natively with
matplotlib — themed via shared.plot_style. Works for both skewer_idealised.html and
skewer_posterior_animation_real.html.

Panels (top→bottom): galaxy strip (colour = inferred class, moving cursor); the three
eigenvalue posterior densities with λ_th marked; the class-probability bar.

Usage:
    render_skewer_video.py <input.html> [--out out.mp4] [--fps 12] [--title "..."]
mp4 (default, ffmpeg) plays inline in Keynote; pass --out ....gif for a looping gif.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from shared.plot_style import apply_style, COSMIC_WEB_COLORS, TEXT_COLOR  # noqa: E402


def load_payload(html_path: Path) -> dict:
    h = html_path.read_text(encoding="utf-8")
    m = re.search(r"var D=(\{.*?\});var F=", h, re.S)
    if not m:
        raise SystemExit(f"Could not find embedded payload (var D=...) in {html_path}")
    return json.loads(m.group(1))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("html", type=Path, help="skewer_*.html to render")
    ap.add_argument("--out", type=Path, default=None, help="output .mp4 (default) or .gif")
    ap.add_argument("--fps", type=int, default=12)
    ap.add_argument("--dpi", type=int, default=110, help="render dpi (lower => smaller gif).")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    d = load_payload(args.html)
    frames = d["frames"]
    xg = np.asarray(d["xg"])
    th = float(d["th"])
    order = d["order"]
    ec = d["ec"]
    gal = d["gal"]
    out = args.out or args.html.with_suffix(".mp4")
    title = args.title or ("Idealised line-of-sight skewer — discrete → continuous"
                           if "idealis" in args.html.name else
                           "Line-of-sight skewer — real DESI NPE posteriors")

    apply_style()
    ymax = max(1e-3, max(max(c) for f in frames for c in f["kde"])) * 1.06
    lname = [r"$\lambda_1$", r"$\lambda_2$", r"$\lambda_3$"]

    fig = plt.figure(figsize=(9, 7.6))
    gs = GridSpec(3, 1, height_ratios=[1.0, 2.6, 0.5], hspace=0.42, figure=fig)
    ax_gal, ax_post, ax_bar = (fig.add_subplot(gs[i]) for i in range(3))
    fig.suptitle(title, fontsize=15, y=0.975)

    # --- galaxy strip (static scatter + moving cursor) ---
    gs_s, gs_t = np.asarray(gal["s"]), np.asarray(gal["t"])
    gcol = [COSMIC_WEB_COLORS[c] for c in gal["c"]]
    ax_gal.scatter(gs_s, gs_t, s=4, c=gcol, alpha=0.55, linewidths=0)
    ax_gal.set_xlim(d["smin"], d["smax"]); ax_gal.set_ylim(d["tmin"], d["tmax"])
    ax_gal.set_xlabel("distance along skewer  [Mpc]"); ax_gal.set_ylabel("transverse [Mpc]")
    ax_gal.set_title("galaxies along the skewer  (colour = inferred class)", fontsize=12)
    cursor = ax_gal.axvline(frames[0]["s"], color=TEXT_COLOR, lw=1.6)

    def update(i):
        f = frames[i]
        cursor.set_xdata([f["s"], f["s"]])

        # posteriors panel
        ax_post.clear()
        ax_post.axvspan(th, xg[-1], color=TEXT_COLOR, alpha=0.05)
        ax_post.axvline(th, ls="--", lw=1.4, color="#9a9a93")
        ax_post.text(th + 0.03, ymax * 0.93, r"$\lambda_{th}=0.2$", color=TEXT_COLOR, fontsize=11)
        for k in range(3):
            c = np.asarray(f["kde"][k])
            ax_post.fill_between(xg, c, color=ec[k], alpha=0.14)
            ax_post.plot(xg, c, color=ec[k], lw=2.1, label=lname[k])
        ax_post.set_xlim(xg[0], xg[-1]); ax_post.set_ylim(0, ymax)
        ax_post.set_xlabel(r"eigenvalue $\lambda$"); ax_post.set_ylabel("posterior density")
        ax_post.legend(loc="upper right", ncol=3, fontsize=11)
        ax_post.text(0.012, 0.95, f"s = {f['s']:.0f} Mpc", transform=ax_post.transAxes,
                     va="top", fontsize=11, color=TEXT_COLOR)

        # class-probability bar
        ax_bar.clear()
        left = 0.0
        for cls in order:
            p = max(0.0, float(f["probs"][cls]))
            if p > 0:
                ax_bar.barh(0, p, left=left, height=1.0, color=COSMIC_WEB_COLORS[cls],
                            edgecolor="#000000", linewidth=0.4)
                if p > 0.07:
                    ax_bar.text(left + p / 2, 0, f"{cls[0].upper()}\n{p*100:.0f}%",
                                ha="center", va="center", fontsize=9, color="#000000")
            left += p
        ax_bar.set_xlim(0, 1); ax_bar.set_ylim(-0.5, 0.5)
        ax_bar.set_yticks([]); ax_bar.set_xticks([])
        for s in ax_bar.spines.values():
            s.set_visible(False)
        ax_bar.set_xlabel("inferred class probability", fontsize=11)
        return ()

    ani = FuncAnimation(fig, update, frames=len(frames), interval=1000 / args.fps, blit=False)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.suffix.lower() == ".gif":
        ani.save(str(out), writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    else:
        # Needs an H.264-capable ffmpeg; falls back to a sibling .gif if encoding fails
        # (e.g. the NERSC system ffmpeg has no libx264).
        try:
            ani.save(str(out), writer=FFMpegWriter(fps=args.fps, bitrate=3000), dpi=args.dpi)
        except Exception as e:
            out = out.with_suffix(".gif")
            print(f"[skewer-video] mp4 encode failed ({type(e).__name__}); writing {out.name} "
                  "(Keynote plays gifs natively).", flush=True)
            ani.save(str(out), writer=PillowWriter(fps=args.fps), dpi=args.dpi)
    plt.close(fig)
    print(f"[skewer-video] wrote {out}  ({len(frames)} frames @ {args.fps} fps)", flush=True)


if __name__ == "__main__":
    main()
