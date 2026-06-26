#!/usr/bin/env python3
"""Mass-anchored figure: massive-halo recovery vs target smoothing, contrasted with the
global R^2(lambda1) optimum. Results from mass_anchored_cluster_test.py (logM>13)."""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

scale = np.array([6, 7, 9, 10, 11, 12, 16, 20])
auc_true = np.array([0.789, 0.770, 0.735, 0.721, 0.708, 0.696, 0.660, 0.637])  # true l1 -> massive halo
auc_pred = np.array([0.722, 0.716, 0.700, 0.694, 0.688, 0.683, 0.671, 0.665])  # model recovers massive
globalR2 = np.array([0.506, 0.547, 0.588, 0.592, 0.588, 0.579, 0.527, 0.466])  # global l1 learnability

apply_style()
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(scale, auc_true, "o-", lw=2.5, color=ACCENT_COLORS["magenta"],
        label=r"AUC: true $\lambda_1$ finds massive halos ($M{>}10^{13}$)")
ax.plot(scale, auc_pred, "s--", lw=2, color="#EB336F",
        label="AUC: model recovers massive halos")
ax.axvline(7, ls=":", color="#9a9a93", label="current 7 Mpc/h")
ax.set_xlabel("target smoothing scale [Mpc/h]")
ax.set_ylabel("massive-halo recovery AUC", color=ACCENT_COLORS["magenta"])
ax.annotate("clusters: finest is best\n(monotonic)", (11, 0.70),
            color=ACCENT_COLORS["magenta"], fontsize=9)

ax2 = ax.twinx()
ax2.plot(scale, globalR2, "^-", lw=2, color=ACCENT_COLORS["blue"], alpha=0.9,
         label=r"global $R^2(\lambda_1)$ (bulk)")
ax2.axvline(10, ls="--", color="#F5C144")
ax2.set_ylabel(r"global $R^2(\lambda_1)$", color=ACCENT_COLORS["blue"])
ax2.annotate("global opt ~10\n(bulk, NOT clusters)", (10.2, 0.50),
             color="#F5C144", fontsize=9)

ax.set_title("Mass-anchored clusters favour FINE smoothing\n(opposite to the global-accuracy optimum)")
l1, lab1 = ax.get_legend_handles_labels(); l2, lab2 = ax2.get_legend_handles_labels()
ax.legend(l1 + l2, lab1 + lab2, loc="lower left", fontsize=8)
fig.tight_layout()
out = Path("/pscratch/sd/d/dkololgi/tng_illustris/figures/smoothing_scale_study")
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "mass_anchored_cluster_recovery.png", bbox_inches="tight", dpi=200)
print("Saved:", out / "mass_anchored_cluster_recovery.png")
