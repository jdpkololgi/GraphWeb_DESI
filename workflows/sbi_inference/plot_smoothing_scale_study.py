#!/usr/bin/env python3
"""Figure for the target-smoothing-scale study (results from
smoothing_scale_investigation.py, wedge footprint, BGS-like downsample).
Left: lambda1/lambda2 learnability + cluster fraction vs smoothing scale.
Right: best-matching density aperture vs smoothing scale (scale-matching).
"""
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
r2_l1 = np.array([0.506, 0.547, 0.588, 0.592, 0.588, 0.579, 0.527, 0.466])
r2_l2 = np.array([0.631, 0.645, 0.638, 0.623, 0.605, 0.583, 0.490, 0.397])
cfrac = np.array([0.085, 0.066, 0.040, 0.032, 0.022, 0.015, 0.003, 0.000])
best_ap = np.array([7, 10, 14, 14, 14, 14, 20, 20])
# cluster-conditioned recovery (the metric that actually matters for clusters)
complete = np.array([0.553, 0.549, 0.544, 0.518, 0.479, 0.470, 0.169, 0.000])

apply_style()
fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 5))

axL.plot(scale, r2_l1, "o-", lw=2, color=ACCENT_COLORS["blue"], label=r"global $R^2(\lambda_1)$ (bulk)")
axL.plot(scale, complete, "o-", lw=2.5, color=ACCENT_COLORS["magenta"], label="cluster completeness")
axL.axvline(7, ls=":", color="#9a9a93", label="current 7 Mpc/h")
axL.axvline(10, ls="--", color="#F5C144", label=r"global-$R^2$ opt ~10")
axL.annotate("global accuracy ↑", (12, 0.58), color=ACCENT_COLORS["blue"], fontsize=8)
axL.annotate("clusters ↓", (12.5, 0.40), color=ACCENT_COLORS["magenta"], fontsize=8)
axL.set_xlabel("target smoothing scale [Mpc/h]"); axL.set_ylabel("score")
axL.set_title("Global accuracy vs CLUSTER recovery (anti-correlated)")
ax2 = axL.twinx()
ax2.plot(scale, cfrac, "^-", lw=1.2, color="#A1FCDD", alpha=0.7)
ax2.set_ylabel(r"cluster fraction $(\lambda_1>0.2)$", color="#A1FCDD")
ax2.tick_params(axis="y", colors="#A1FCDD")
axL.legend(loc="center left", fontsize=8)

# cross-scale fate: clusters defined at finest scale, tracked outward
phys = np.array([0.747, 0.415, 0.327, 0.218, 0.140, 0.023, 0.002])  # rs7..rs20 persistence
sc2 = np.array([7, 9, 10, 11, 12, 16, 20])
axR.plot(sc2, phys, "o-", lw=2, color=ACCENT_COLORS["magenta"], label="still a cluster (physical)")
axR.axhline(0.5, ls=":", color="#9a9a93"); axR.axvline(10, ls="--", color="#F5C144")
axR.annotate("⅔ of fine-scale\nclusters gone by 10", (10.3, 0.55), color="#F2F2F2", fontsize=8)
axR.set_xlabel("target smoothing scale [Mpc/h]")
axR.set_ylabel("frac of rs6 clusters surviving")
axR.set_title("Compact clusters dissolve under smoothing")
axR.legend(fontsize=9)

fig.suptitle("T-web target smoothing scale study (BGS-like sampling, wedge footprint)", y=1.04)
fig.tight_layout()
out = Path("/pscratch/sd/d/dkololgi/tng_illustris/figures/smoothing_scale_study")
out.mkdir(parents=True, exist_ok=True)
fig.savefig(out / "smoothing_scale_learnability.png", bbox_inches="tight", dpi=200)
print("Saved:", out / "smoothing_scale_learnability.png")
