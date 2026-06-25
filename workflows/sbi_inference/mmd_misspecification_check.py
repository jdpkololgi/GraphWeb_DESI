#!/usr/bin/env python3
"""Summary-space MMD model-misspecification check (Schmitt et al. 2023 style),
as referenced by the BayesFlow tutorial ("summary space MMD check for model
misspecification").

In our amortized NPE the GNN encoder IS the summary network, so the 80-d node
embeddings are the summary space. Misspecification on DESI = the DESI embedding
distribution drifting away from the Abacus (training) embedding distribution.

We compute an unbiased MMD^2 with a median-heuristic RBF kernel and a permutation
null. Crucially we also compute the Abacus-vs-Abacus split-half MMD as the
"well-specified" floor: with ~thousands of points any real shift gives p~0, so
the interpretable quantity is how far DESI sits ABOVE that in-distribution floor.

No GPU: operates on already-saved embeddings.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist


def mmd2_unbiased(X, Y, gamma):
    """Unbiased MMD^2 with RBF kernel k(a,b)=exp(-gamma||a-b||^2)."""
    Kxx = np.exp(-gamma * cdist(X, X, "sqeuclidean"))
    Kyy = np.exp(-gamma * cdist(Y, Y, "sqeuclidean"))
    Kxy = np.exp(-gamma * cdist(X, Y, "sqeuclidean"))
    m, n = len(X), len(Y)
    np.fill_diagonal(Kxx, 0.0); np.fill_diagonal(Kyy, 0.0)
    return (Kxx.sum() / (m * (m - 1)) + Kyy.sum() / (n * (n - 1))
            - 2.0 * Kxy.mean())


def median_gamma(Z):
    """Median-heuristic bandwidth on a pooled sample Z -> gamma = 1/(2 sigma^2)."""
    sub = Z[np.random.choice(len(Z), min(2000, len(Z)), replace=False)]
    d = cdist(sub, sub, "euclidean")
    sigma = np.median(d[d > 0])
    return 1.0 / (2.0 * sigma ** 2)


def perm_pvalue(X, Y, gamma, obs, n_perm=200):
    Z = np.vstack([X, Y]); m = len(X)
    null = np.empty(n_perm)
    for i in range(n_perm):
        idx = np.random.permutation(len(Z))
        null[i] = mmd2_unbiased(Z[idx[:m]], Z[idx[m:]], gamma)
    p = (1 + np.sum(null >= obs)) / (1 + n_perm)
    return p, null


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--abacus-npz", type=Path, required=True)
    ap.add_argument("--desi-npz", type=Path, required=True)
    ap.add_argument("--n", type=int, default=2000, help="equal subsample size per group")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    np.random.seed(args.seed)

    ab = np.load(args.abacus_npz)["embeddings"].astype(np.float64)
    de = np.load(args.desi_npz)["embeddings"].astype(np.float64)
    print(f"Abacus embeddings {ab.shape}; DESI embeddings {de.shape}")

    # standardize by ABACUS (training) statistics -> the units the model expects
    mu, sd = ab.mean(0), ab.std(0) + 1e-12
    ab = (ab - mu) / sd; de = (de - mu) / sd

    floor, cross, ratios, pvals = [], [], [], []
    for r in range(args.repeats):
        ia = np.random.choice(len(ab), 2 * args.n, replace=False)
        A1, A2 = ab[ia[:args.n]], ab[ia[args.n:]]          # disjoint Abacus halves
        D = de[np.random.choice(len(de), args.n, replace=False)]
        gamma = median_gamma(np.vstack([A1, D]))
        f = mmd2_unbiased(A1, A2, gamma)                    # well-specified floor
        c = mmd2_unbiased(A1, D, gamma)                     # Abacus vs DESI
        p, _ = perm_pvalue(A1, D, gamma, c, n_perm=args.n_perm)
        floor.append(f); cross.append(c); ratios.append(c / max(f, 1e-12)); pvals.append(p)
        print(f"  repeat {r}: floor(Ab|Ab)={f:.5f}  cross(Ab|DESI)={c:.5f}  "
              f"ratio={c/max(f,1e-12):.1f}x  p={p:.4f}")

    floor, cross = np.asarray(floor), np.asarray(cross)
    # Effect size: how far the DESI cross-MMD^2 sits above the in-distribution
    # floor, in units of the floor's own scatter (the floor mean is ~0 for an
    # unbiased estimator, so a ratio is meaningless; a z-score is not).
    floor_scale = max(floor.std(), 1e-9)
    z = (cross.mean() - floor.mean()) / floor_scale
    mmd_dist = np.sqrt(max(cross.mean(), 0.0))   # actual MMD (not squared), standardized 80-d space
    print("\n=== summary-space MMD misspecification check ===")
    print(f"  well-specified floor  MMD^2(Abacus,Abacus) = {floor.mean():+.5f} +/- {floor.std():.5f}")
    print(f"  observed              MMD^2(Abacus,DESI)   = {cross.mean():+.5f} +/- {cross.std():.5f}")
    print(f"  MMD distance (sqrt, standardized 80-d)     = {mmd_dist:.4f}")
    print(f"  effect size (DESI above floor)             = {z:.1f} sigma_floor")
    print(f"  permutation p-value (median)               = {np.median(pvals):.4f}")
    detected = np.median(pvals) < 0.05
    severe = mmd_dist > 0.2   # ~half a full distribution swap on standardized embeddings
    verdict = ("SEVERE misspecification: large summary-space shift" if (detected and severe)
               else "DETECTABLE but small shift: significant yet modest magnitude (typical sim->real)"
               if detected else
               "WELL-SPECIFIED: DESI summaries indistinguishable from the Abacus floor")
    print(f"  verdict: {verdict}")


if __name__ == "__main__":
    main()
