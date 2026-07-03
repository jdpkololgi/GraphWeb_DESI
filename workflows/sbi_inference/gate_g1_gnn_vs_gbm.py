#!/usr/bin/env python3
"""Gate G1 — GNN-vs-GBM capacity check (roadmap v2, Track 1).

Question: does the trained Attentional GraphNetwork extract more eigenvalue
information than a gradient-boosted tree given the IDENTICAL node features and
splits? GNN >> GBM => capacity/message-passing matters; GNN ~= GBM => the
hand-crafted feature set is the binding constraint (information-limited).

No new training: the GNN side is the SI production model's posterior-mean
predictions on the test split (abacus_self_flowjax_preds.npz). The GBM side is
HistGBM trained on the cache's train split (graph.nodes = the exact scaled
features the GNN consumed) and scored on the same test rows.

GO criterion (adopt capacity work) : GNN lambda1 R^2 - GBM lambda1 R^2 > ~0.03.
CPU only.
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import r2_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, required=True, help="SI training cache pkl")
    ap.add_argument("--self-eval-npz", type=Path, required=True,
                    help="abacus_self_flowjax_preds.npz (GNN posterior means on test split)")
    ap.add_argument("--lambda-th", type=float, default=0.2)
    args = ap.parse_args()

    cache = pickle.load(open(args.cache, "rb"))
    X = np.asarray(cache["graph"].nodes, np.float64)          # model-input features
    eig = np.asarray(cache["eigenvalues_raw"], np.float64)     # [N,3] raw truth
    train, val, test = (np.asarray(m).astype(bool) for m in cache["masks"])
    d = np.load(args.self_eval_npz)
    node_index = d["node_index"].astype(int)
    gnn_mean = d["lambda_mean"].astype(np.float64)
    gnn_truth = d["eig_truth"].astype(np.float64)

    # alignment guards: npz test rows must be the cache's test mask, same truth
    ti = np.where(test)[0]
    assert np.array_equal(np.sort(node_index), np.sort(ti)), "test-row mismatch"
    order = np.argsort(node_index)
    ni = node_index[order]; gm = gnn_mean[order]; gt = gnn_truth[order]
    assert np.allclose(gt, eig[ni], atol=1e-4), "truth mismatch between npz and cache"
    print(f"aligned: {len(ni)} test rows; features {X.shape[1]}-d; train={train.sum()}")

    print(f"\n{'':10s}  {'GBM R2':>8s}  {'GNN R2':>8s}  {'delta':>7s}   (test split, raw eigenvalues)")
    deltas = {}
    for k, nm in enumerate(["lambda1", "lambda2", "lambda3"]):
        gbm = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.08, max_depth=6)
        gbm.fit(X[train], eig[train, k])
        pg = gbm.predict(X[ni])
        r_gbm = r2_score(eig[ni, k], pg)
        r_gnn = r2_score(eig[ni, k], gm[:, k])
        deltas[nm] = r_gnn - r_gbm
        print(f"{nm:10s}  {r_gbm:8.3f}  {r_gnn:8.3f}  {r_gnn - r_gbm:+7.3f}")
        if k == 0:
            # cluster-conditioned slice (true lambda1 > lambda_th)
            clu = eig[ni, 0] > args.lambda_th
            if clu.sum() > 50:
                print(f"{'  cluster':10s}  {r2_score(eig[ni,0][clu], pg[clu]):8.3f}"
                      f"  {r2_score(eig[ni,0][clu], gm[clu,0]):8.3f}"
                      f"   (n={clu.sum()}; Spearman GBM {spearmanr(eig[ni,0][clu], pg[clu]).statistic:+.2f}"
                      f" / GNN {spearmanr(eig[ni,0][clu], gm[clu,0]).statistic:+.2f})")

    go = deltas["lambda1"] > 0.03
    print(f"\nGATE G1: delta(lambda1) = {deltas['lambda1']:+.3f} -> "
          f"{'GO (capacity matters; message passing adds signal)' if go else 'NO-GO (information-limited at this feature set)'}")


if __name__ == "__main__":
    main()
