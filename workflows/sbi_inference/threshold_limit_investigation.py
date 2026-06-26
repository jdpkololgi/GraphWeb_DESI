#!/usr/bin/env python3
"""Why is lambda1 unpredictable near lambda_th? (Phase-0 variant D, eigenvalue space.)

Tests the candidate causes of the boundary-zone limit found in variant C:
  [1] FEATURE-SCALE MISMATCH: the geom features are Delaunay-local (~few Mpc) but the
      target lambda1 is the 7 Mpc/h-smoothed tidal field. Does adding a SCALE-MATCHED
      aperture density (count within 7 & 10 Mpc/h) lift lambda1 R^2? (the lever might be
      density-at-the-right-scale, not kinematics).
  [2] IRREDUCIBLE FLOOR: for galaxies with near-identical features, how much does true
      lambda1 still scatter? median local std(lambda1) / global std = the noise floor
      the features cannot beat. Computed for geom7 and for the scale-matched set.
  [3] METRIC ARTIFACT vs GENUINE: in the threshold zone report RMSE & Spearman (not just
      R^2, which is degenerate in a narrow-variance slice).
  [4] SOFT LABEL / lambda_th sensitivity: how true-cluster fraction and a geom-regressor's
      cluster<->filament confusion vary with lambda_th -> motivates reporting P(lambda1>lambda_th).

Mock truth only. CPU.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table
from scipy.spatial import cKDTree
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score

from velocity_dispersion_precheck import los_dispersion, fog_aniso
from velocity_dispersion_aperture_precheck import aperture_dispersion


def cvpred(X, y, folds):
    m = HistGradientBoostingRegressor(max_iter=250, learning_rate=0.08, max_depth=6)
    return cross_val_predict(m, X, y, cv=folds)


def floor(Xstd, y, k=30, sample=8000, seed=0):
    """median local std(y) among k nearest neighbours in standardized feature space."""
    rng = np.random.default_rng(seed)
    tree = cKDTree(Xstd)
    s = rng.permutation(len(Xstd))[:sample]
    _, nn = tree.query(Xstd[s], k=k)
    return float(np.median(np.std(y[nn], axis=1)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock-points-xyz", type=Path, required=True)
    ap.add_argument("--mock-gnn-arrays", type=Path, required=True)
    ap.add_argument("--mock-targets-fits", type=Path, required=True)
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-n", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed); h = cosmo.h

    pos = np.load(args.mock_points_xyz).astype(np.float64)
    mx = np.load(args.mock_gnn_arrays)["x"].astype(np.float64)
    t = Table.read(args.mock_targets_fits)
    l1 = np.asarray(t["LAMBDA1"], float)

    # scale-matched aperture features: density (neighbour count) @7,@10 Mpc/h + FoGAniso@10
    feats = {}
    for R in (7.0, 10.0):
        sp, spp, cnt = aperture_dispersion(pos, R / h)
        feats[f"logdens{int(R)}"] = np.log1p(cnt)
        if R == 10.0:
            feats["fog10"] = fog_aniso(sp, spp)
    dens7 = feats["logdens7"]; dens10 = feats["logdens10"]; fog10 = feats["fog10"]

    idx = np.arange(len(pos))
    if len(idx) > args.max_n:
        idx = rng.permutation(idx)[: args.max_n]
    y = l1[idx]; Xg = mx[idx]
    Xd = np.column_stack([Xg, dens7[idx], dens10[idx]])           # +scale-matched density
    Xa = np.column_stack([Xg, dens7[idx], dens10[idx], fog10[idx]])  # +density +kinematics
    glob_std = y.std()
    print(f"mock n={len(idx)}; global std(lambda1)={glob_std:.3f}; lambda_th={args.lambda_th}")

    # [1] scale-matched density vs kinematics for lambda1
    print("\n[1] feature-scale-match test -- CV R^2(lambda1)")
    for tag, X in [("geom7 (local)", Xg), ("geom7+dens@7,10 (scale-matched)", Xd),
                   ("geom7+dens+FoG@10", Xa)]:
        print(f"    {tag:34s} R^2 = {r2_score(y, cvpred(X, y, args.folds)):.3f}")

    # [2] irreducible floor (median local std lambda1 / global)
    print("\n[2] irreducible floor -- median local std(lambda1)/global (1.0=unpredictable, 0=determined)")
    for tag, X in [("geom7", Xg), ("geom7+dens+FoG", Xa)]:
        Xs = (X - X.mean(0)) / (X.std(0) + 1e-9)
        f = floor(Xs, y, k=30, sample=8000, seed=args.seed)
        print(f"    {tag:16s} floor={f:.3f}  ratio={f/glob_std:.2f}")

    # [3] threshold-zone metric (RMSE & Spearman, not R^2)
    print("\n[3] threshold zone (|lambda1-lambda_th|<0.15) -- RMSE & Spearman vs R^2")
    pg = cvpred(Xg, y, args.folds); pa = cvpred(Xa, y, args.folds)
    z = np.abs(y - args.lambda_th) < 0.15
    for tag, p in [("geom7", pg), ("geom7+dens+FoG", pa)]:
        rmse = np.sqrt(np.mean((y[z] - p[z]) ** 2)); sp = spearmanr(y[z], p[z]).statistic
        print(f"    {tag:16s} zone(n={z.sum()}) RMSE={rmse:.3f}  Spearman={sp:+.3f}  "
              f"(zone std(lambda1)={y[z].std():.3f})")

    # [4] lambda_th sensitivity: true cluster fraction + geom-regressor confusion
    print("\n[4] lambda_th sensitivity (cluster = lambda1>lambda_th; geom regressor)")
    print("    lambda_th  true_cluster_frac  pred_cluster_frac  cluster->non(miss)  Spearman(all)")
    sp_all = spearmanr(y, pg).statistic
    for th in (0.0, 0.1, 0.2, 0.3, 0.4):
        tc = y > th; pc = pg > th
        miss = float(np.mean(~pc[tc])) if tc.any() else float("nan")
        print(f"    {th:6.1f}     {tc.mean():.3f}              {pc.mean():.3f}             "
              f"{miss:.3f}            {sp_all:+.3f}")

    print("\nReadout: if [1] dens@scale lifts R^2 a lot -> the limit was FEATURE SCALE, not "
          "kinematics. If [2] ratio stays ~1 and [3] Spearman~0 -> genuinely irreducible near "
          "threshold. [4] shows why soft P(lambda1>lambda_th) is the honest product.")


if __name__ == "__main__":
    main()
