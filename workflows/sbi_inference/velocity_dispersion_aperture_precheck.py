#!/usr/bin/env python3
"""Phase-0 variant B: LOS velocity dispersion at FIXED APERTURES matched to the
T-web target smoothing scale (7 Mpc/h).

The original precheck measured dispersion over Delaunay neighbours (~few Mpc) and
found no cluster signal. But the targets (lambda1/2/3) are the tidal eigenvalues of
the density field Gaussian-smoothed at 7 Mpc/h, so the matched kinematic scale is
~7 Mpc/h, not the virial/Delaunay scale. This sweeps aperture radius (in Mpc/h)
around 7 and asks, on MOCK TRUTH (CWEB), whether the dispersion features then
separate true clusters from filaments and add cluster-recovery headroom.

Aperture neighbours come from scipy cKDTree.query_pairs(r) -> fed to the same
los_dispersion() bincount accumulator. CPU only.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import roc_auc_score, recall_score, f1_score

from velocity_dispersion_precheck import los_dispersion, fog_aniso  # reuse
from astropy.table import Table


def aperture_dispersion(pos, radius):
    """sigma_par, sigma_perp, count over all neighbours within `radius` (cKDTree pairs)."""
    tree = cKDTree(pos)
    pairs = tree.query_pairs(radius, output_type="ndarray")  # [P,2], each undirected once
    if len(pairs) == 0:
        z = np.zeros(len(pos)); return z, z, z.astype(int)
    return los_dispersion(pos, pairs.T)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock-points-xyz", type=Path, required=True)
    ap.add_argument("--mock-gnn-arrays", type=Path, required=True, help="for 7 geom features (x)")
    ap.add_argument("--mock-targets-fits", type=Path, required=True)
    ap.add_argument("--radii-hmpc", type=float, nargs="+", default=[3, 5, 7, 10, 14])
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-n", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)
    h = cosmo.h
    print(f"Planck18 h={h:.4f}; apertures (Mpc/h -> Mpc): " +
          ", ".join(f"{r:g}->{r/h:.1f}" for r in args.radii_hmpc))

    pos = np.load(args.mock_points_xyz).astype(np.float64)
    mx = np.load(args.mock_gnn_arrays)["x"].astype(np.float64)
    cweb = np.asarray(Table.read(args.mock_targets_fits)["CWEB"]).astype(int)
    print(f"mock nodes {len(pos)}; CWEB counts {np.bincount(cweb)}")

    midx = np.arange(len(pos))
    if len(midx) > args.max_n:
        midx = rng.permutation(midx)[: args.max_n]
    yc = cweb[midx]

    def cv4(X, tag):
        clf = HistGradientBoostingClassifier(max_iter=250, learning_rate=0.08, max_depth=6,
                                             class_weight="balanced")
        pred = cross_val_predict(clf, X, yc, cv=args.folds)
        rec = recall_score(yc, pred, labels=[0, 1, 2, 3], average=None, zero_division=0)
        f1m = f1_score(yc, pred, average="macro", zero_division=0)
        c2f = float(np.mean(pred[yc == 3] == 2)) if (yc == 3).any() else float("nan")
        print(f"      {tag:18s} cluster_recall={rec[3]:.3f}  filament_recall={rec[2]:.3f}  "
              f"macroF1={f1m:.3f}  cluster->filament={c2f:.3f}")
        return rec[3], c2f

    print("\n=== geom-only baseline (class-balanced HistGBM) ===")
    cr0, c2f0 = cv4(mx[midx], "geom-only(7)")

    for R in args.radii_hmpc:
        r_mpc = R / h
        spar, sperp, cnt = aperture_dispersion(pos, r_mpc)
        fog = fog_aniso(spar, sperp)
        delta = spar**2 - sperp**2          # signed LOS anisotropy (FoG>0 / Kaiser<0)
        ok = cnt >= 5
        cols = {"sigma_par": spar, "sigma_perp": sperp, "FoGAniso": fog, "delta": delta}
        print(f"\n=== aperture {R:g} Mpc/h ({r_mpc:.1f} Mpc); median neighbours={np.median(cnt):.0f} ===")
        # univariate true cluster(3) vs filament(2)
        mm = ok & np.isin(cweb, [2, 3])
        yt = (cweb[mm] == 3).astype(int)
        for nm in ("sigma_par", "FoGAniso", "delta"):
            ft = cols[nm][mm]; a = roc_auc_score(yt, ft); a = max(a, 1 - a)
            print(f"      [uni] {nm:10s} AUC={a:.3f}  "
                  f"(clu med {np.median(ft[yt==1]):+.3f} vs fil med {np.median(ft[yt==0]):+.3f})")
        # incremental cluster recovery for candidate feature subsets
        subsets = {"delta": ["delta"], "FoGAniso": ["FoGAniso"],
                   "FoG+delta": ["FoGAniso", "delta"],
                   "all4": ["sigma_par", "sigma_perp", "FoGAniso", "delta"]}
        for tag, names in subsets.items():
            feat = np.stack([cols[c] for c in names], axis=1)[midx]
            cr1, c2f1 = cv4(np.hstack([mx[midx], feat]), f"geom+{tag}")
            print(f"        Δrecall {cr1-cr0:+.3f}   Δ(cluster->filament) {c2f1-c2f0:+.3f}")

    print("\nGATE: a matched aperture passes if cluster recall rises / cluster->filament "
          "drops materially vs geom-only, and a univariate AUC clearly exceeds 0.5.")


if __name__ == "__main__":
    main()
