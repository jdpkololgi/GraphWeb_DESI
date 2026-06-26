#!/usr/bin/env python3
"""Phase-0 variant C: does the scale-matched FoGAniso feature improve prediction of
the EIGENVALUES (the actual regression targets), especially lambda1?

The NPE regresses (lambda1,lambda2,lambda3); the T-web class is just a threshold
(lambda_th=0.2) applied afterward. cluster vs filament == whether lambda1 crosses
lambda_th. So the right pre-check is in eigenvalue space, not on discrete labels:

  [A] incremental CV R^2 for lambda1/2/3 from geom7 vs geom7+FoGAniso(aperture).
  [B] does FoGAniso explain the geometry's lambda1 RESIDUAL (the info geometry
      misses)? corr + 1-D R^2 of FoGAniso -> (lambda1_true - lambda1_pred_geom).
  [C] same incremental R^2 but restricted to the THRESHOLD ZONE (lambda1 near 0.2),
      where resolving lambda1 actually flips the class.

Mock truth LAMBDA1/2/3 are the 7 Mpc/h-smoothed tidal eigenvalues. CPU only.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score
from scipy.stats import pearsonr

from velocity_dispersion_precheck import los_dispersion, fog_aniso
from velocity_dispersion_aperture_precheck import aperture_dispersion


def cv_r2(X, y, folds):
    m = HistGradientBoostingRegressor(max_iter=250, learning_rate=0.08, max_depth=6)
    return r2_score(y, cross_val_predict(m, X, y, cv=folds)), cross_val_predict(m, X, y, cv=folds)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock-points-xyz", type=Path, required=True)
    ap.add_argument("--mock-gnn-arrays", type=Path, required=True)
    ap.add_argument("--mock-targets-fits", type=Path, required=True)
    ap.add_argument("--radii-hmpc", type=float, nargs="+", default=[8, 10, 12])
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-n", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    h = cosmo.h

    pos = np.load(args.mock_points_xyz).astype(np.float64)
    mx = np.load(args.mock_gnn_arrays)["x"].astype(np.float64)
    t = Table.read(args.mock_targets_fits)
    lam = np.stack([np.asarray(t[f"LAMBDA{i}"], float) for i in (1, 2, 3)], axis=1)  # [N,3]
    print(f"mock nodes {len(pos)}; lambda1 range [{lam[:,0].min():.2f},{lam[:,0].max():.2f}] "
          f"median {np.median(lam[:,0]):.3f}; lambda_th={args.lambda_th}")

    idx = np.arange(len(pos))
    if len(idx) > args.max_n:
        idx = rng.permutation(idx)[: args.max_n]
    Xg = mx[idx]
    Y = lam[idx]

    # geometry-only baselines (once)
    print("\n=== geometry-only CV R^2 (baseline) ===")
    base_r2, base_pred = {}, {}
    for j, nm in enumerate(["lambda1", "lambda2", "lambda3"]):
        r2, pred = cv_r2(Xg, Y[:, j], args.folds)
        base_r2[nm] = r2; base_pred[nm] = pred
        print(f"    geom -> {nm}: R^2 = {r2:.3f}")
    resid_l1 = Y[:, 0] - base_pred["lambda1"]   # geometry's lambda1 error

    for R in args.radii_hmpc:
        r_mpc = R / h
        spar, sperp, cnt = aperture_dispersion(pos, r_mpc)
        fog = fog_aniso(spar, sperp)[idx]
        print(f"\n=== aperture {R:g} Mpc/h ({r_mpc:.1f} Mpc) ===")
        # [A] incremental R^2 per eigenvalue
        for j, nm in enumerate(["lambda1", "lambda2", "lambda3"]):
            r2f, _ = cv_r2(np.column_stack([Xg, fog]), Y[:, j], args.folds)
            print(f"    [A] {nm}: geom R^2 {base_r2[nm]:.3f} -> geom+FoG {r2f:.3f}  (Δ {r2f-base_r2[nm]:+.3f})")
        # [B] does FoG explain geometry's lambda1 residual?
        pr = pearsonr(fog, resid_l1)[0]
        r2_resid, _ = cv_r2(fog.reshape(-1, 1), resid_l1, args.folds)
        print(f"    [B] FoG -> lambda1 residual: pearson r={pr:+.3f}, 1-D R^2={r2_resid:.3f} "
              f"({'carries missing info' if r2_resid > 0.02 else 'negligible'})")
        # [C] threshold zone: lambda1 within +/-0.15 of lambda_th
        zone = np.abs(Y[:, 0] - args.lambda_th) < 0.15
        if zone.sum() > 1000:
            rg, _ = cv_r2(Xg[zone], Y[zone, 0], args.folds)
            rf, _ = cv_r2(np.column_stack([Xg, fog])[zone], Y[zone, 0], args.folds)
            print(f"    [C] lambda1 in threshold zone (n={zone.sum()}): "
                  f"geom R^2 {rg:.3f} -> geom+FoG {rf:.3f}  (Δ {rf-rg:+.3f})")

    print("\nGATE: FoG helps if [A] Δ R^2(lambda1) clearly >0, [B] explains the residual, "
          "and [C] improves lambda1 specifically in the threshold zone.")


if __name__ == "__main__":
    main()
