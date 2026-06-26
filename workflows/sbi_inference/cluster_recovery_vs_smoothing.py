#!/usr/bin/env python3
"""Does the global R^2(lambda1) optimum at ~10 Mpc/h actually mean BETTER CLUSTER
inference? Physics critique (JDPK): clusters are compact, so larger smoothing washes
them out -> global accuracy can rise while cluster recovery falls. Global R^2 (3-6%
clusters) is dominated by the non-cluster bulk and is the wrong objective for clusters.

Cluster-CONDITIONED metrics vs target smoothing scale (cutsky, BGS-like downsample):
  [A] AUC(true-cluster vs rest) of predicted lambda1  -- threshold/shrinkage-robust
      separability of clusters. If this peaks at SMALLER scale than global R^2 -> the
      10 Mpc/h optimum is misleading for clusters.
  [B] completeness@true-rate: of the top-N predicted lambda1 (N=#true clusters), frac
      that are true clusters -- rank-based recall decoupled from shrinkage.
  [C] cluster-regime bias: mean(pred-true) lambda1 for true clusters (shrinkage).
  [D] cross-scale fate: clusters defined at the FINEST scale -> what frac are still
      clusters (physical persistence) and still recoverable (model) at larger smoothing.

Reuses aperture-density features (fixed) from smoothing_scale_investigation. CPU.
"""
from __future__ import annotations
import argparse
import numpy as np
import fitsio
from astropy.cosmology import Planck18 as cosmo
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score, roc_auc_score

from smoothing_scale_investigation import SCALE_FILES, radec_z_to_xyz, aperture_density


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", nargs=2, type=float, default=[120, 160])
    ap.add_argument("--dec", nargs=2, type=float, default=[14.5, 30.6])
    ap.add_argument("--zr", nargs=2, type=float, default=[0.2, 0.3])
    ap.add_argument("--target-n", type=int, default=120000)
    ap.add_argument("--apertures-hmpc", type=float, nargs="+", default=[3, 5, 7, 10, 14, 20, 28])
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed); h = cosmo.h; th = args.lambda_th

    base = fitsio.read(SCALE_FILES[7], columns=["RA", "DEC", "Z"])
    ra, dec, z = base["RA"], base["DEC"], base["Z"]
    fp = ((ra >= args.ra[0]) & (ra < args.ra[1]) & (dec >= args.dec[0]) & (dec < args.dec[1]) &
          (z >= args.zr[0]) & (z < args.zr[1]))
    fidx = np.where(fp)[0]
    if len(fidx) > args.target_n:
        fidx = fidx[rng.permutation(len(fidx))[: args.target_n]]
    pos = radec_z_to_xyz(ra[fidx], dec[fidx], z[fidx])
    dens = np.column_stack([np.log1p(aperture_density(pos, pos, R / h)) for R in args.apertures_hmpc])
    print(f"galaxies={len(fidx)}; lambda_th={th}\n")

    def cvpred(y):
        m = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, max_depth=6)
        return cross_val_predict(m, dens, y, cv=args.folds)

    scales = sorted(SCALE_FILES)
    l1_by_s, pred_by_s = {}, {}
    print("scale  cfrac  globalR2  | AUC(cluster)  complete@rate  cluster_bias")
    for s in scales:
        l1 = fitsio.read(SCALE_FILES[s], columns=["LAMBDA1"])["LAMBDA1"][fidx].astype(np.float64)
        pred = cvpred(l1)
        l1_by_s[s] = l1; pred_by_s[s] = pred
        clu = l1 > th; nclu = int(clu.sum())
        gR2 = r2_score(l1, pred)
        if nclu > 20:
            auc = roc_auc_score(clu, pred)
            topN = np.argsort(pred)[::-1][:nclu]          # top-N predicted (N=#true clusters)
            complete = float(clu[topN].mean())            # rank-based recall=precision@N
            bias = float(np.mean(pred[clu] - l1[clu]))    # shrinkage in cluster regime
            print(f"rs{s:<3d}  {clu.mean():.3f}  {gR2:.3f}    | {auc:.3f}        {complete:.3f}"
                  f"          {bias:+.3f}")
        else:
            print(f"rs{s:<3d}  {clu.mean():.3f}  {gR2:.3f}    | (too few clusters: {nclu})")

    # [D] cross-scale fate: clusters at the finest scale, tracked outward
    ref = scales[0]
    C = l1_by_s[ref] > th
    print(f"\n[D] cross-scale fate of clusters defined at rs{ref} (N={int(C.sum())}):")
    print("    scale  still_cluster(phys)  recovered_by_model(top-N@ref-rate)")
    nref = int(C.sum())
    for s in scales:
        phys = float((l1_by_s[s][C] > th).mean())                  # physical persistence
        topN = np.argsort(pred_by_s[s])[::-1][:nref]               # model's top-N predictions
        rec = np.zeros(len(C), bool); rec[topN] = True
        recovered = float(rec[C].mean())                            # frac of ref-clusters in model top-N
        print(f"    rs{s:<3d}   {phys:.3f}                {recovered:.3f}")

    print("\nReadout: if AUC(cluster)/completeness PEAK at a SMALLER scale than the global-R^2 "
          "peak (~10), then optimising global R^2 hurts clusters (JDPK's concern confirmed). "
          "[D] quantifies how fast true clusters dissolve & become unrecoverable with smoothing.")


if __name__ == "__main__":
    main()
