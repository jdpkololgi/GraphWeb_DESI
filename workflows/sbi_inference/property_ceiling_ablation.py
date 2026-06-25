#!/usr/bin/env python3
"""Property-information ceiling for the spatial-only NPE (option 2 scoping).

Question: are we limiting the model by using features completely determined by
the spatial distribution of galaxies? A clean, non-circular bound is REDUNDANCY:
how well does the geometry the model ALREADY uses (its 80-d GNN embedding) predict
each FastSpecFit galaxy property?

  geometry predicts property well  -> property is redundant with geometry
                                      -> adding it as a model input buys little.
  geometry predicts property poorly-> property carries INDEPENDENT information
                                      -> real headroom from adding it (but you'd
                                         first have to paint that property onto
                                         the Abacus mocks to train on it).

Why this and not "do properties improve environment recovery": DESI has no TRUE
environment label, so that direct test is impossible here; predicting the model's
own (geometry-derived) class from properties is circular. Redundancy is the
honest, available proxy.

CPU only. Uses the closure-test parquet (properties) + preds npz (embeddings).
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.model_selection import cross_val_predict, cross_val_score
from sklearn.metrics import r2_score, roc_auc_score

CONT = ["LOGMSTAR", "gr", "log_sSFR", "DN4000"]
BIN = ["quenched"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds-npz", type=Path, required=True)
    ap.add_argument("--closure-parquet", type=Path, required=True)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-n", type=int, default=40000, help="subsample for speed")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    d = np.load(args.preds_npz)
    emb = pd.DataFrame(d["embeddings"])
    emb.columns = [f"e{i}" for i in range(emb.shape[1])]
    emb["global_node_id"] = d["global_node_id"]
    emb["hc_pred"] = d["hard_class"]

    prop = pd.read_parquet(args.closure_parquet)
    df = emb.merge(prop, on="global_node_id", how="inner")
    df = df[df["has_props"]].copy()
    print(f"joined galaxies with FastSpecFit props: {len(df)}")
    if len(df) > args.max_n:
        df = df.iloc[rng.permutation(len(df))[: args.max_n]].copy()
        print(f"subsampled to {len(df)} for speed")

    E = [c for c in df.columns if c.startswith("e") and c[1:].isdigit()]
    X = df[E].to_numpy(np.float64)

    def reg_redundancy(mask=None, tag="all"):
        sub = df if mask is None else df[mask]
        Xs = sub[E].to_numpy(np.float64)
        print(f"\n--- regression redundancy ({tag}, n={len(sub)}) ---")
        print("  property     geometry->property R^2   reading")
        for col in CONT:
            y = sub[col].to_numpy(np.float64)
            ok = np.isfinite(y)
            if ok.sum() < 500:
                print(f"  {col:10s}  (too few finite values)"); continue
            m = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, max_depth=6)
            yp = cross_val_predict(m, Xs[ok], y[ok], cv=args.folds)
            r2 = r2_score(y[ok], yp)
            verdict = "redundant (low headroom)" if r2 > 0.5 else \
                      "partly independent" if r2 > 0.2 else "mostly INDEPENDENT (headroom)"
            print(f"  {col:10s}  {r2:6.3f}                  {verdict}")

    def clf_redundancy(mask=None, tag="all"):
        sub = df if mask is None else df[mask]
        Xs = sub[E].to_numpy(np.float64)
        print(f"\n--- classification redundancy ({tag}, n={len(sub)}) ---")
        for col in BIN:
            y = sub[col].to_numpy(int)
            ok = np.isfinite(sub[col].to_numpy(float)) & np.isin(y, [0, 1])
            if ok.sum() < 500 or len(np.unique(y[ok])) < 2:
                print(f"  {col:10s}  (insufficient)"); continue
            m = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.08, max_depth=6)
            pp = cross_val_predict(m, Xs[ok], y[ok], cv=args.folds, method="predict_proba")[:, 1]
            auc = roc_auc_score(y[ok], pp)
            base = y[ok].mean()
            verdict = "redundant (low headroom)" if auc > 0.8 else \
                      "partly independent" if auc > 0.65 else "mostly INDEPENDENT (headroom)"
            print(f"  {col:10s}  geometry->{col} AUC={auc:.3f} (base rate {base:.2f})  {verdict}")

    print("\n================ PROPERTY-INFORMATION CEILING ================")
    print("How much of each galaxy property does the model's geometry already know?")
    reg_redundancy(tag="all galaxies")
    clf_redundancy(tag="all galaxies")

    # cluster-focused: the failure mode lives in dense/cluster-classified galaxies
    clu = df["hc_pred"].to_numpy() == 3
    if clu.sum() > 1000:
        reg_redundancy(mask=clu, tag="model-cluster galaxies")
        clf_redundancy(mask=clu, tag="model-cluster galaxies")

    print("\nNOTE: high R^2/AUC = geometry already encodes the property (low ceiling).")
    print("      low  R^2/AUC = property is independent info the spatial model lacks")
    print("      (= real headroom, but requires painting it onto the Abacus mocks to train).")


if __name__ == "__main__":
    main()
