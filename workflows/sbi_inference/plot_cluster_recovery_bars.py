#!/usr/bin/env python3
"""Slide-9 'this didn't come for free' bar chart: DESI cluster-fraction recovery
toward Abacus truth across the domain-correction sequence.

Baseline → inference-time edge-adapt → node+edge adapt → SI retrain (production),
with the Abacus-truth cluster fraction as a reference line. Numbers are read live
from each run's summary.json (no hardcoding); the % of the gap closed by the
production SI model is annotated.

Themed via shared.plot_style; written into the SI run dir by default.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from shared.plot_style import apply_style, finalize_axes, COSMIC_WEB_COLORS  # noqa: E402
from shared.config_paths import GRAPHWEB_SCRATCH_ROOT  # noqa: E402

# (run-dir suffix, short label) in recovery order.
STAGES = [
    ("desi_wedge_flowjax_linear", "Baseline"),
    ("desi_wedge_flowjax_linear_edgeadapt", "Edge-adapt"),
    ("desi_wedge_flowjax_linear_fulladapt", "Node+edge\nadapt"),
    ("desi_wedge_flowjax_linear_si", "Per-graph norm.\n(production)"),
]


def _cluster_frac(summary_path: Path) -> float:
    s = json.loads(summary_path.read_text())
    return float(s["desi"]["class_fractions_npe"]["cluster"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    flx = Path(f"{GRAPHWEB_SCRATCH_ROOT}/flowjax_inference_outputs")
    ap.add_argument("--inference-root", type=Path, default=flx)
    ap.add_argument("--outdir", type=Path,
                    default=flx / "desi_wedge_flowjax_linear_si",
                    help="where to write cluster_recovery_bars.png (default: the SI run dir).")
    args = ap.parse_args()

    apply_style()
    vals, labels = [], []
    for suffix, label in STAGES:
        sp = args.inference_root / suffix / "summary.json"
        if not sp.exists():
            print(f"[recovery] skip (no summary): {sp}", flush=True)
            continue
        vals.append(_cluster_frac(sp)); labels.append(label)
    vals = np.array(vals)
    # Abacus-truth cluster fraction (reference) from the SI run's reference_fractions.
    si_summary = json.loads((args.inference_root / "desi_wedge_flowjax_linear_si" / "summary.json").read_text())
    truth = float(si_summary["reference_fractions"]["abacus_truth"]["cluster"])

    gold = COSMIC_WEB_COLORS["cluster"]
    alphas = np.linspace(0.5, 1.0, len(vals))
    fig, ax = plt.subplots(figsize=(9, 6))
    x = np.arange(len(vals))
    for xi, v, a in zip(x, vals, alphas):
        ax.bar(xi, v, 0.62, color=gold, alpha=a, edgecolor="#F2F2F2", lw=0.8)
        ax.text(xi, v + 0.0015, f"{v:.3f}", ha="center", va="bottom", fontsize=11)
    ax.axhline(truth, ls="--", lw=1.6, color="#F2F2F2")
    ax.text(len(vals) - 0.5, truth + 0.0015, f"Abacus truth = {truth:.3f}",
            ha="right", va="bottom", fontsize=11, color="#F2F2F2")

    gap_closed = (vals[-1] - vals[0]) / (truth - vals[0]) * 100.0 if len(vals) >= 2 else float("nan")
    ax.annotate(f"{gap_closed:.0f}% of the gap closed",
                xy=(x[-1], vals[-1]), xytext=(x[-1] - 1.1, truth - 0.004),
                fontsize=12, color=gold,
                arrowprops=dict(arrowstyle="->", color=gold, lw=1.4))

    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_ylim(0.0, max(truth, vals.max()) * 1.18)
    finalize_axes(ax, r"DESI cluster-fraction recovery toward truth ($\lambda_{th}=0.2$)",
                  "domain-correction stage", "inferred cluster fraction", legend=False)
    out = args.outdir / "cluster_recovery_bars.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight"); plt.close(fig)
    print(f"[recovery] wrote {out}", flush=True)
    print(f"[recovery] sequence: {dict(zip([l.replace(chr(10),' ') for l in labels], np.round(vals,3)))}; "
          f"truth={truth:.3f}; gap closed {gap_closed:.0f}%", flush=True)


if __name__ == "__main__":
    main()
