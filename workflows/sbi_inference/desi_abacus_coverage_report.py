"""Training-support coverage report in the SI model-input space.

Quantifies how far DESI (and DESI cluster-candidate) galaxies fall beyond the
Abacus training support, per node feature, in the EXACT space the SI model
ingests: per-graph-median normalise scale cols [0,2,3,4,5,6] (Clustering col 1
excluded), then the cache box-cox PowerTransformer. Beyond-support fraction =
the extrapolation the flow is asked to do (drives clusters -> filaments).

CPU only.
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path
import numpy as np

NAMES = ["Degree", "Clustering", "Density", "NeighDensity", "I_eig1", "I_eig2", "I_eig3"]
SI_COLS = [0, 2, 3, 4, 5, 6]   # scale-carrying cols normalised per-graph (matches SI cache/inference)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--si-cache", type=Path, required=True)
    ap.add_argument("--desi-gnn-arrays", type=Path, required=True)
    ap.add_argument("--preds-npz", type=Path, required=True, help="for DESI hard_class (cluster candidates)")
    args = ap.parse_args()

    cache = pickle.load(open(args.si_cache, "rb"))
    nfs = cache["node_feature_scaler"]
    A = np.asarray(cache["graph"].nodes, np.float64)          # Abacus model-input space (SI+boxcox)
    train = np.asarray(cache["masks"][0]).astype(bool)
    A_tr = A[train]                                            # training support

    x_raw = np.load(args.desi_gnn_arrays)["x"].astype(np.float64).copy()
    for c in SI_COLS:                                          # per-graph-median SI normalise (DESI's own median)
        x_raw[:, c] /= float(np.median(x_raw[:, c]))
    D = nfs.transform(x_raw + 1e-6)                            # into the literal model-input space

    hard = np.load(args.preds_npz)["hard_class"]
    clu = hard == 3 if len(hard) == len(D) else np.zeros(len(D), bool)

    print(f"Abacus train nodes {A_tr.shape[0]}; DESI nodes {D.shape[0]}; "
          f"DESI model-cluster {int(clu.sum())}")
    print("\nfeature        train_p99.9   frac_DESI>p99.9   frac_DESI>train_max   cluster-cand frac>p99.9")
    for k in range(7):
        tp = np.percentile(A_tr[:, k], 99.9); tmax = A_tr[:, k].max()
        f99 = np.mean(D[:, k] > tp); fmax = np.mean(D[:, k] > tmax)
        fclu = np.mean(D[clu, k] > tp) if clu.any() else np.nan
        flag = "  <== OOD" if (f99 > 0.02 or fmax > 0.005) else ""
        print(f"{NAMES[k]:13s}  {tp:8.3f}    {f99:9.4f}        {fmax:9.4f}           {fclu:9.4f}{flag}")

    print("\nReading: large frac>train_max in the density-family cols (Density/NeighDensity)")
    print("for cluster candidates = the flow extrapolating beyond training support there.")


if __name__ == "__main__":
    main()
