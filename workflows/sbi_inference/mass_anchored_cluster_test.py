#!/usr/bin/env python3
"""Mass-anchored cluster test (definitive): anchor 'cluster' to the SCALE-INDEPENDENT
physical label HALO_MASS (a direct column in the master cutsky catalog; no fragile
CompaSO index join), and ask:
  [sanity] does HALO_MASS correlate with lambda1 and with the galaxy density field?
  [A] which TARGET smoothing scale's true lambda1 best identifies massive halos?
      (AUC of true lambda1(s) for M>M_th) -- the natural cluster scale, mass-anchored.
  [B] does the model (aperture-density features) recover massive-halo galaxies, and
      how does that depend on smoothing? AUC/completeness of predicted lambda1(s).
  [C] features -> mass directly: is mass recoverable from the galaxy distribution?
  [D] alignment of T-web cluster (lambda1(s)>lambda_th) with mass cluster vs s.

Master cutsky is row-aligned with the rs6-24 lambda catalogs (same 63.9M galaxies).
"""
from __future__ import annotations
import argparse
import numpy as np
import fitsio
from astropy.cosmology import Planck18 as cosmo
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score, r2_score

from smoothing_scale_investigation import SCALE_FILES, radec_z_to_xyz, aperture_density

MASTER = "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_23032026/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs.fits"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", nargs=2, type=float, default=[120, 160])
    ap.add_argument("--dec", nargs=2, type=float, default=[14.5, 30.6])
    ap.add_argument("--zr", nargs=2, type=float, default=[0.2, 0.3])
    ap.add_argument("--target-n", type=int, default=120000)
    ap.add_argument("--apertures-hmpc", type=float, nargs="+", default=[3, 7, 10, 14])
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--mass-th-log", type=float, nargs="+", default=[12.5, 13.0, 13.5],
                    help="log10(Msun/h) cluster thresholds (HALO_MASS column is in 1e10 Msun/h)")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed); h = cosmo.h; th = args.lambda_th

    m = fitsio.read(MASTER, columns=["RA", "DEC", "Z", "HALO_MASS", "FILE_NUM"])
    ra, dec, z, M = m["RA"], m["DEC"], m["Z"], m["HALO_MASS"].astype(np.float64)
    fp = ((ra >= args.ra[0]) & (ra < args.ra[1]) & (dec >= args.dec[0]) & (dec < args.dec[1]) &
          (z >= args.zr[0]) & (z < args.zr[1]))
    fidx = np.where(fp)[0]
    if len(fidx) > args.target_n:
        fidx = fidx[rng.permutation(len(fidx))[: args.target_n]]
    pos = radec_z_to_xyz(ra[fidx], dec[fidx], z[fidx])
    Mf = M[fidx]
    # HALO_MASS is in 1e10 Msun/h -> log10(Msun/h) = log10(Mf) + 10
    logM = np.log10(np.where(Mf > 0, Mf, np.nan)) + 10.0
    print(f"footprint galaxies={len(fidx)}; HALO_MASS(1e10 Msun/h) median={np.nanmedian(Mf):.2f}; "
          f"log10(M/[Msun/h]) median={np.nanmedian(logM):.2f}, "
          f"range [{np.nanmin(logM):.1f},{np.nanmax(logM):.1f}], frac_valid={np.mean(np.isfinite(logM)):.3f}")

    dens = np.column_stack([np.log1p(aperture_density(pos, pos, R / h)) for R in args.apertures_hmpc])

    # row-aligned lambda1 at each smoothing
    scales = sorted(SCALE_FILES)
    l1 = {s: fitsio.read(SCALE_FILES[s], columns=["LAMBDA1"])["LAMBDA1"][fidx].astype(np.float64) for s in scales}

    # SANITY GATE
    print("\n[sanity] (correct join => positive)")
    print(f"  Spearman(logM, lambda1@rs7) = {spearmanr(logM, l1[7], nan_policy='omit').statistic:+.3f}")
    yM = (logM > 13.0).astype(int)
    aucM7 = roc_auc_score(yM, cross_val_predict(
        HistGradientBoostingClassifier(max_iter=200, learning_rate=0.08, max_depth=6),
        dens, yM, cv=args.folds, method="predict_proba")[:, 1])
    print(f"  features -> P(logM>13) AUC = {aucM7:.3f}  ({'OK' if aucM7 > 0.6 else 'SUSPECT'})")

    for mthlog in args.mass_th_log:
        mclu = logM > mthlog; nM = int(mclu.sum())
        print(f"\n=== mass cluster logM>{mthlog:.1f} (M>{10**mthlog:.1e} Msun/h): frac={mclu.mean():.4f} (N={nM}) ===")
        print("  scale  cfrac(l1>th)  AUC_true_l1(mass)  AUC_pred_l1(mass)  complete@rate  Tweb∩mass/Tweb")
        for s in scales:
            ls = l1[s]; tweb = ls > th
            auc_true = roc_auc_score(mclu, ls)
            pred = cross_val_predict(HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, max_depth=6),
                                     dens, ls, cv=args.folds)
            auc_pred = roc_auc_score(mclu, pred)
            topN = np.argsort(pred)[::-1][:nM]; comp = float(mclu[topN].mean())
            purity = float((tweb & mclu).sum() / max(tweb.sum(), 1))
            print(f"  rs{s:<3d}  {tweb.mean():.3f}        {auc_true:.3f}             {auc_pred:.3f}"
                  f"            {comp:.3f}         {purity:.3f}")

    print("\nReadout: [A] AUC_true_l1(mass) peak = smoothing where the T-web label best matches "
          "real massive halos. [B] AUC_pred/completeness = can the model recover them, vs scale. "
          "[D] Tweb∩mass/Tweb = purity of the T-web cluster label as a massive-halo finder.")


if __name__ == "__main__":
    main()
