#!/usr/bin/env python3
"""Gates G1.5 (RSD penalty) + G2 (luminosity weighting) — roadmap v2, Track 1.

Master cutsky carries BOTH Z (observed, RSD) and Z_COSMO (cosmological) plus
R_MAG_ABS and rs7 LAMBDA truth (sampled at the REAL-SPACE host-halo x_com —
verified in annotate_cutsky_with_tweb_eigs.py, so the targets carry no RSD
sampling offset and the G1.5 penalty is purely input-side information loss).

For ONE downsampled BGS-like sample in the wedge footprint we build aperture
features from a 2x2 grid: positions {z-space, real-space} x weighting
{count, +luminosity}, and measure GBM CV R^2(lambda1) + mass-anchored cluster
recovery for each cell.

  G1.5 RSD penalty   = R2(real, count) - R2(z, count)   [upper bound on ALL
                       LOS-aware / equivariant machinery, this feature family]
  G2   lum gain      = R2(z, count+lum) - R2(z, count)  [GO if > ~0.03]

CPU only.
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

MASTER = ("/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_23032026/"
          "cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs.fits")


def radec_z_to_xyz(ra, dec, z):
    zt = np.linspace(0, max(float(np.max(z)), 0.7), 4000)
    dt = cosmo.comoving_distance(zt).value
    d = np.interp(z, zt, dt)
    r, dd = np.deg2rad(ra), np.deg2rad(dec)
    return np.vstack([d*np.cos(dd)*np.cos(r), d*np.cos(dd)*np.sin(r), d*np.sin(dd)]).T


def aperture_features(pos, radii_mpc, lum=None):
    """log(1+count) per aperture; if lum given, also log(1e-9+sum L) per aperture."""
    tree = cKDTree(pos)
    cols = []
    for r in radii_mpc:
        nb = tree.query_ball_point(pos, r, workers=-1)
        cnt = np.fromiter((len(ix) for ix in nb), dtype=np.int64, count=len(pos))
        cols.append(np.log1p(cnt))
        if lum is not None:
            ls = np.fromiter((lum[ix].sum() for ix in nb), dtype=np.float64, count=len(pos))
            cols.append(np.log(1e-9 + ls))
    return np.column_stack(cols)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", nargs=2, type=float, default=[120, 160])
    ap.add_argument("--dec", nargs=2, type=float, default=[14.5, 30.6])
    ap.add_argument("--zr", nargs=2, type=float, default=[0.2, 0.3])
    ap.add_argument("--target-n", type=int, default=120000)
    ap.add_argument("--apertures-hmpc", type=float, nargs="+", default=[3, 7, 10, 14])
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--mass-log-th", type=float, default=13.0)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed); h = cosmo.h

    m = fitsio.read(MASTER, columns=["RA", "DEC", "Z", "Z_COSMO", "R_MAG_ABS",
                                     "HALO_MASS", "LAMBDA1"])
    fp = ((m["RA"] >= args.ra[0]) & (m["RA"] < args.ra[1]) &
          (m["DEC"] >= args.dec[0]) & (m["DEC"] < args.dec[1]) &
          (m["Z"] >= args.zr[0]) & (m["Z"] < args.zr[1]) &
          np.isfinite(m["LAMBDA1"]))
    fidx = np.where(fp)[0]
    if len(fidx) > args.target_n:
        fidx = fidx[rng.permutation(len(fidx))[: args.target_n]]
    sub = m[fidx]
    l1 = sub["LAMBDA1"].astype(np.float64)
    logM = np.log10(np.maximum(sub["HALO_MASS"].astype(np.float64), 1e-3)) + 10.0
    lum = 10.0 ** (-0.4 * sub["R_MAG_ABS"].astype(np.float64))   # relative L_r
    pos_z = radec_z_to_xyz(sub["RA"], sub["DEC"], sub["Z"])
    pos_r = radec_z_to_xyz(sub["RA"], sub["DEC"], sub["Z_COSMO"])
    rsd_disp = np.linalg.norm(pos_z - pos_r, axis=1)
    print(f"sample n={len(sub)}; median |RSD displacement| = {np.median(rsd_disp):.2f} Mpc; "
          f"frac(logM>{args.mass_log_th:g}) = {(logM > args.mass_log_th).mean():.4f}")

    radii = [R / h for R in args.apertures_hmpc]
    mclu = logM > args.mass_log_th; nM = int(mclu.sum())
    tclu = l1 > args.lambda_th

    def score(X, tag):
        gbm = HistGradientBoostingRegressor(max_iter=250, learning_rate=0.08, max_depth=6)
        pred = cross_val_predict(gbm, X, l1, cv=args.folds)
        r2 = r2_score(l1, pred)
        auc = roc_auc_score(mclu, pred) if nM > 50 else float("nan")
        topN = np.argsort(pred)[::-1][:nM]
        comp = float(mclu[topN].mean()) if nM > 50 else float("nan")
        zone = np.abs(l1 - args.lambda_th) < 0.15
        r2z = r2_score(l1[zone], pred[zone]) if zone.sum() > 500 else float("nan")
        print(f"  {tag:24s} R2(l1)={r2:.3f}  massAUC={auc:.3f}  mass-compl@rate={comp:.3f}  zoneR2={r2z:+.3f}")
        return r2

    print("\n=== 2x2 grid: positions x weighting (GBM CV, lambda1) ===")
    r_zc = score(aperture_features(pos_z, radii), "z-space, count")
    r_zl = score(aperture_features(pos_z, radii, lum=lum), "z-space, count+lum")
    r_rc = score(aperture_features(pos_r, radii), "real-space, count")
    r_rl = score(aperture_features(pos_r, radii, lum=lum), "real-space, count+lum")

    pen = r_rc - r_zc
    gain = r_zl - r_zc
    print(f"\nGATE G1.5: RSD penalty (count feats)     = {pen:+.3f}  -> "
          f"{'LARGE: LOS-aware/equivariant machinery has headroom (G4 motivated)' if pen > 0.05 else 'SMALL: little for RSD-aware architecture to recover'}")
    print(f"GATE G2  : luminosity gain (z-space)     = {gain:+.3f}  -> "
          f"{'GO: include flux-weighted features' if gain > 0.03 else 'NO-GO: luminosity adds little'}")
    print(f"           interaction (real+lum vs real) = {r_rl - r_rc:+.3f}")
    print(f"           context: true-cluster frac (l1>th) = {tclu.mean():.3f}")


if __name__ == "__main__":
    main()
