#!/usr/bin/env python3
"""Why is the environment<->property signal weak? Decompose physics vs inference noise.

Tests:
 1. Independent local density (kNN in redshift-space comoving coords, n(z)-controlled)
    -> does the INFERRED trace track real local density? (model quality)
 2. property vs {inferred trace, independent density, mass} Spearman -> is the
    property-density relation itself weak (physics) or only the inferred one (noise)?
 3. Model accuracy on the mock: predicted vs TRUE trace (abacus self-eval).
 4. Uncertainty domination: per-galaxy posterior width vs galaxy-to-galaxy spread.
 5. Mass-controlled correlations + top/bottom-decile effect sizes.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from astropy.cosmology import Planck18

CIG = Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
           "desi_wedge_cigale_hz/desi_wedge_env_props.parquet")
ASELF = Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
             "abacus_self_linear_si/abacus_self_flowjax_preds.npz")


def sr(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return spearmanr(a[m], b[m]).correlation, int(m.sum())


def main():
    df = pd.read_parquet(CIG)
    n = len(df)
    # ---- 1. independent local density (redshift-space comoving xyz) ----
    chi = Planck18.comoving_distance(df["z"].to_numpy()).value  # Mpc
    ra, dec = np.deg2rad(df["ra"].to_numpy()), np.deg2rad(df["dec"].to_numpy())
    xyz = np.c_[chi*np.cos(dec)*np.cos(ra), chi*np.cos(dec)*np.sin(ra), chi*np.sin(dec)]
    tree = cKDTree(xyz)
    zbin = pd.qcut(df["z"], 25, labels=False, duplicates="drop")
    for k in (5, 10):
        dk, _ = tree.query(xyz, k=k+1, workers=-1)
        logdens = np.log10(k / (4/3*np.pi*dk[:, -1]**3))
        med = pd.Series(logdens).groupby(zbin).transform("median").to_numpy()
        df[f"ddens{k}"] = logdens - med          # n(z)-controlled log density contrast
    print(f"[data] {n:,} wedge galaxies; kNN density built (n(z)-controlled)")

    ms = df["has_ssfr"].to_numpy()
    tr = df["trace_lambda"].to_numpy()
    dd = df["ddens10"].to_numpy()
    print("\n=== KEY: does the INFERRED trace track independent local density? ===")
    r, _ = sr(tr, dd); print(f"  Spearman(inferred trace, kNN density contrast) = {r:+.3f}")
    r5, _ = sr(tr, df['ddens5'].to_numpy()); print(f"  (k=5 density) = {r5:+.3f}")

    print("\n=== property vs INFERRED trace  vs  INDEPENDENT density  vs  mass ===")
    props = {"log_sSFR": df["log_sSFR"].to_numpy(),
             "gr": df["gr"].to_numpy(),
             "quenched": df["quenched"].to_numpy().astype(float),
             "LOGMSTAR": df["LOGMSTAR"].to_numpy()}
    lm = df["LOGMSTAR"].to_numpy()
    print(f"  {'property':10s} {'vs trace':>10s} {'vs kNN dens':>12s} {'vs mass':>9s}")
    for name, v in props.items():
        v2 = np.where(ms, v, np.nan) if name != "LOGMSTAR" else v
        rt, _ = sr(v2, tr); rd, _ = sr(v2, dd); rm, _ = sr(v2, lm)
        print(f"  {name:10s} {rt:>+10.3f} {rd:>+12.3f} {rm:>+9.3f}")

    print("\n=== mass-controlled (within logM* tertiles): property vs trace / density ===")
    mt = pd.qcut(pd.Series(lm), 3, labels=False)
    for name in ("log_sSFR", "quenched"):
        v = props[name]
        row = []
        for t in range(3):
            sel = ms & (mt.to_numpy() == t)
            rt, _ = sr(v[sel], tr[sel]); rd, _ = sr(v[sel], dd[sel])
            row.append(f"T{t}: tr {rt:+.3f}/dens {rd:+.3f}")
        print(f"  {name:10s} " + "  ".join(row))

    print("\n=== effect size: top vs bottom DECILE (Δ median property) ===")
    for xn, xv in (("inferred trace", tr), ("kNN density", dd)):
        lo, hi = np.nanpercentile(xv[ms], [10, 90])
        for name in ("log_sSFR", "gr"):
            v = props[name]
            dmed = np.nanmedian(v[ms & (xv > hi)]) - np.nanmedian(v[ms & (xv < lo)])
            print(f"  {xn:14s} Δmedian {name:8s} (top10%-bot10%) = {dmed:+.3f}")

    print("\n=== 3. model accuracy on the MOCK (predicted vs TRUE trace, in-dist) ===")
    a = np.load(ASELF)
    tt = a["eig_truth"].sum(1); tp = a["lambda_mean"].sum(1)
    R2 = 1 - np.sum((tt-tp)**2)/np.sum((tt-tt.mean())**2)
    rho = spearmanr(tt, tp).correlation
    print(f"  trace: R^2 = {R2:.3f}  Spearman = {rho:.3f}  (attenuation factor ~ {rho:.2f})")
    for i in range(3):
        tti, tpi = a["eig_truth"][:, i], a["lambda_mean"][:, i]
        print(f"  lambda{i+1}: R^2 = {1-np.sum((tti-tpi)**2)/np.sum((tti-tti.mean())**2):.3f}"
              f"  Spearman = {spearmanr(tti, tpi).correlation:.3f}")

    print("\n=== 4. uncertainty domination on DESI (posterior width vs spread) ===")
    for i in (1, 2, 3):
        spread = np.std(df[f"lambda_{i}_mean"])
        width = np.median(df[f"lambda_{i}_std"])
        print(f"  lambda{i}: galaxy-spread(mean)={spread:.3f}  median posterior width={width:.3f}"
              f"  ratio width/spread={width/spread:.2f}")
    tr_spread = np.std(tr)
    tr_width = np.median(np.sqrt(df["lambda_1_std"]**2 + df["lambda_2_std"]**2 + df["lambda_3_std"]**2))
    print(f"  trace: spread(mean)={tr_spread:.3f}  ~post width={tr_width:.3f}  ratio={tr_width/tr_spread:.2f}")

    # ---- diagnostic figure ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    sys.path.insert(0, "/global/homes/d/dkololgi/TNG/Illustris")
    from shared.plot_style import apply_style, ACCENT_COLORS
    apply_style()
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2), layout="constrained")
    hb = ax[0].hexbin(dd[ms], tr[ms], gridsize=55, mincnt=3, cmap="viridis", bins="log")
    ax[0].set_xlabel("independent kNN density contrast"); ax[0].set_ylabel(r"inferred trace $E[\Sigma\lambda]$")
    ax[0].set_title(f"Model recovers real density  (ρ = {sr(tr, dd)[0]:+.2f})", fontsize=13)
    ax[0].set_xlim(np.nanpercentile(dd[ms], [1, 99])); ax[0].set_ylim(np.nanpercentile(tr[ms], [1, 99]))
    labels = [r"log sSFR", r"$(g-r)$", "quenched"]
    keys = ["log_sSFR", "gr", "quenched"]
    rm = [abs(sr(props[k] if k != "quenched" else props[k], lm)[0]) for k in keys]
    rd = [abs(sr(np.where(ms, props[k], np.nan), dd)[0]) for k in keys]
    rtr = [abs(sr(np.where(ms, props[k], np.nan), tr)[0]) for k in keys]
    yy = np.arange(3); h = 0.26
    ax[1].barh(yy+h, rm, h, color="#F2F2F2", label="vs stellar mass")
    ax[1].barh(yy, rd, h, color=ACCENT_COLORS["blue"], label="vs independent density")
    ax[1].barh(yy-h, rtr, h, color=ACCENT_COLORS["magenta"], label="vs inferred trace")
    ax[1].set_yticks(yy); ax[1].set_yticklabels(labels); ax[1].set_xlabel(r"$|\rho_{\rm Spearman}|$")
    ax[1].set_title("Mass dominates; environment weak but real\n(inferred ≈ independent density)", fontsize=13)
    ax[1].legend(fontsize=10, loc="lower right", frameon=False)
    out = CIG.parent.parent.parent / "figures/desi_wedge_cigale_hz/env_signal_diagnosis.png"
    fig.savefig(out, dpi=200, bbox_inches="tight"); fig.savefig(str(out).replace(".png", ".pdf"), bbox_inches="tight")
    print("\nSaved figure:", out)


if __name__ == "__main__":
    main()
