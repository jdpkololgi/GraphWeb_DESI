#!/usr/bin/env python3
"""Phase-0 gate for the LOS velocity-dispersion node feature (plan: keen-twirling-prism).

Computes, per galaxy, the line-of-sight decomposition of graph-neighbour
displacements (reusing the bincount pattern from plot_fog_los_alignment.py):

  sigma_par  = sqrt(<(d . r_hat)^2>)              # spread along the LOS
  sigma_perp = sqrt(0.5 * <|d - (d.r_hat) r_hat|^2>)  # transverse spread
  FoGAniso   = sigma_par / sigma_perp             # >1 => collapsed/virialised (FoG)

Then it answers the three gating questions, BEFORE any retrain:
  (2) signal:       does FoGAniso separate model-cluster vs model-filament galaxies?
  (3) redundancy:   does the geometry the model already uses (80-d embedding) predict
                    FoGAniso? low R^2 => the feature is new information.
  (4) transfer:     do the mock (path1) and DESI FoGAniso distributions match? a large
                    gap flags the redshift-error injection branch of the plan.

No new observable: uses the redshift-space comoving positions already built
(observer at origin => r_hat = pos/|pos|). CPU only.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
from astropy.cosmology import Planck18 as cosmo
from astropy.table import Table
from sklearn.ensemble import HistGradientBoostingRegressor, HistGradientBoostingClassifier
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score, roc_auc_score, recall_score, f1_score


def radec_z_to_xyz(ra, dec, z):
    d = cosmo.comoving_distance(np.asarray(z)).value  # Mpc
    r, dd = np.deg2rad(ra), np.deg2rad(dec)
    return np.vstack([d * np.cos(dd) * np.cos(r),
                      d * np.cos(dd) * np.sin(r),
                      d * np.sin(dd)]).T


def los_dispersion(pos, edge_index):
    """Return (sigma_par, sigma_perp, count) per node from neighbour displacements.

    Each undirected edge contributes to BOTH endpoints, each using its own LOS
    direction r_hat = pos/|pos|. Mirrors los_alignment()'s bincount accumulation.
    """
    N = len(pos)
    s, r = edge_index[0], edge_index[1]
    los = pos / (np.linalg.norm(pos, axis=1, keepdims=True) + 1e-12)
    d = pos[r] - pos[s]                       # [E,3] displacement (s -> r)
    d2 = np.einsum("ij,ij->i", d, d)          # |d|^2
    # endpoint s: project d onto r_hat_s ; endpoint r: project (-d) onto r_hat_r
    proj_s = np.einsum("ij,ij->i", d, los[s])
    proj_r = np.einsum("ij,ij->i", d, los[r])  # sign irrelevant (squared)
    par = np.bincount(s, proj_s**2, N) + np.bincount(r, proj_r**2, N)
    perp = (np.bincount(s, d2 - proj_s**2, N) + np.bincount(r, d2 - proj_r**2, N))
    cnt = np.bincount(s, minlength=N) + np.bincount(r, minlength=N)
    safe = np.maximum(cnt, 1)
    sigma_par = np.sqrt(np.maximum(par / safe, 0.0))
    sigma_perp = np.sqrt(np.maximum(0.5 * perp / safe, 0.0))
    return sigma_par, sigma_perp, cnt


def fog_aniso(sigma_par, sigma_perp):
    return sigma_par / np.maximum(sigma_perp, 1e-9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--desi-preds-npz", type=Path, required=True)
    ap.add_argument("--desi-gnn-arrays", type=Path, required=True)
    ap.add_argument("--mock-points-xyz", type=Path, required=True)
    ap.add_argument("--mock-gnn-arrays", type=Path, required=True)
    ap.add_argument("--mock-targets-fits", type=Path, required=True, help="mock CWEB truth, aligned to points_xyz")
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--max-n", type=int, default=40000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed)

    # ---- DESI ----
    d = np.load(args.desi_preds_npz)
    pos = radec_z_to_xyz(d["ra"], d["dec"], d["z"])
    ei = np.load(args.desi_gnn_arrays)["edge_index"]
    spar, sperp, cnt = los_dispersion(pos, ei)
    fog = fog_aniso(spar, sperp)
    ok = cnt >= 3
    hard = d["hard_class"]; emb = d["embeddings"]
    print(f"DESI nodes {len(pos)} (deg>=3: {ok.sum()}); "
          f"FoGAniso median {np.median(fog[ok]):.3f}  sigma_par median {np.median(spar[ok]):.3f} Mpc")

    # (3) redundancy: geometry (embedding) -> FoGAniso ; low R^2 => new info  (DESI, valid)
    print("\n[3] redundancy -- geometry(80-d embedding) -> FoGAniso (low R^2 = NEW info)")
    idx = np.where(ok)[0]
    if len(idx) > args.max_n:
        idx = rng.permutation(idx)[: args.max_n]
    reg = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, max_depth=6)
    yp = cross_val_predict(reg, emb[idx], fog[idx], cv=args.folds)
    r2 = r2_score(fog[idx], yp)
    print(f"    embedding -> FoGAniso  R^2 = {r2:.3f}  "
          f"({'NEW info (headroom)' if r2 < 0.3 else 'partly redundant' if r2 < 0.6 else 'redundant'})")

    # ---- mock (path1): TRUTH available -> the decisive tests ----
    mpos = np.load(args.mock_points_xyz).astype(np.float64)
    mx = np.load(args.mock_gnn_arrays)["x"].astype(np.float64)          # 7 geometric features
    mei = np.load(args.mock_gnn_arrays)["edge_index"]
    cweb = np.asarray(Table.read(args.mock_targets_fits)["CWEB"]).astype(int)  # 0/1/2/3 truth
    mspar, msperp, mcnt = los_dispersion(mpos, mei)
    mfog = fog_aniso(mspar, msperp)
    mok = mcnt >= 3
    veld = np.stack([mspar, msperp, mfog], axis=1)                       # 3 velocity features

    # [2] DECISIVE signal: do velocity features separate TRUE cluster(3) vs filament(2)?
    print("\n[2] TRUE-label signal (mock) -- velocity features, cluster(3) vs filament(2)")
    mm = mok & np.isin(cweb, [2, 3])
    yt = (cweb[mm] == 3).astype(int)
    for name, feat in [("sigma_par", mspar[mm]), ("sigma_perp", msperp[mm]), ("FoGAniso", mfog[mm])]:
        a = roc_auc_score(yt, feat); a = max(a, 1 - a)
        print(f"    {name:10s} AUC={a:.3f}  (cluster med {np.median(feat[yt==1]):.3f} vs "
              f"filament med {np.median(feat[yt==0]):.3f})")

    # [2b] DECISIVE incremental value: cheap 4-class proxy, geom-only vs geom+velocity
    print("\n[2b] incremental cluster recovery (mock truth, HistGBM proxy for the GNN)")
    midx = np.where(mok)[0]
    if len(midx) > args.max_n:
        midx = rng.permutation(midx)[: args.max_n]
    yc = cweb[midx]
    def cv4(X, tag):
        clf = HistGradientBoostingClassifier(max_iter=250, learning_rate=0.08, max_depth=6)
        pred = cross_val_predict(clf, X, yc, cv=args.folds)
        rec = recall_score(yc, pred, labels=[0, 1, 2, 3], average=None, zero_division=0)
        f1m = f1_score(yc, pred, average="macro", zero_division=0)
        # cluster->filament confusion: of true clusters, frac predicted filament
        tc = yc == 3
        c2f = float(np.mean(pred[tc] == 2)) if tc.any() else float("nan")
        print(f"    {tag:16s} cluster_recall={rec[3]:.3f}  filament_recall={rec[2]:.3f}  "
              f"macroF1={f1m:.3f}  cluster->filament={c2f:.3f}")
        return rec[3], c2f
    cr0, c2f0 = cv4(mx[midx], "geom-only(7)")
    cr1, c2f1 = cv4(np.hstack([mx[midx], veld[midx]]), "geom+veld(10)")
    print(f"    --> cluster recall {cr0:.3f} -> {cr1:.3f} (delta {cr1-cr0:+.3f}); "
          f"cluster->filament {c2f0:.3f} -> {c2f1:.3f} (delta {c2f1-c2f0:+.3f})")

    # (4) transfer: mock vs DESI FoGAniso distribution
    print("\n[4] transfer -- mock(path1) vs DESI FoGAniso distribution (large gap => inject z-errors)")
    qs = [10, 25, 50, 75, 90, 99]
    dq = np.percentile(fog[ok], qs); mq = np.percentile(mfog[mok], qs)
    print("    pctl:   " + "  ".join(f"{q:>3d}" for q in qs))
    print("    DESI:   " + "  ".join(f"{v:.2f}" for v in dq))
    print("    mock:   " + "  ".join(f"{v:.2f}" for v in mq))
    med_shift = (np.median(fog[ok]) - np.median(mfog[mok])) / np.median(mfog[mok])
    print(f"    median fractional shift (DESI-mock)/mock = {med_shift:+.1%}  "
          f"({'SMALL (SI ok)' if abs(med_shift) < 0.1 else 'LARGE -> consider z-error injection'})")

    print("\nGATE: proceed to Phase 1 only if [2b] shows a real cluster-recall gain / "
          "cluster->filament drop on mock truth, AND [3] confirms the feature is new info.")


if __name__ == "__main__":
    main()
