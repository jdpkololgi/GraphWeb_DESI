#!/usr/bin/env python3
"""How does the T-web TARGET smoothing scale affect lambda1/lambda2 (and their
learnability from the galaxy distribution)?

Uses the cutsky BGS multi-smoothing eigenvalue catalogs (same 63.9M galaxies, rs
6-24 Mpc/h). Features (galaxy geometry) are FIXED; only the target smoothing varies.
In a wedge-sized footprint, downsampled to BGS-like density, we measure per smoothing
scale s:
  - distribution: std(lambda1), frac(lambda1>lambda_th)
  - predictability: aperture-density features -> lambda1/lambda2 CV R^2
  - scale-matching: which single aperture density best predicts lambda1(s)?
    (hypothesis: aperture ~ s)

fitsio reads only the needed columns. CPU; high-memory node.
"""
from __future__ import annotations
import argparse
import numpy as np
import fitsio
from astropy.cosmology import Planck18 as cosmo
from scipy.spatial import cKDTree
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_val_predict
from sklearn.metrics import r2_score

# scale (Mpc/h) -> catalog path (cutsky, _15d set; same galaxies across scales)
SCALE_FILES = {
    6:  "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_20260427_rsmooth_6/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs6_ngrid2048_thr0p2_15d.fits",
    7:  "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_20260427_rsmooth_7/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs7_ngrid2048_thr0p2_15d.fits",
    9:  "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_20260427_rsmooth_9/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs9_ngrid2048_thr0p2_15d.fits",
    10: "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_20260427_rsmooth_10/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs10_ngrid2048_thr0p2_15d.fits",
    11: "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_20260427_rsmooth_11/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs11_ngrid2048_thr0p2_15d.fits",
    12: "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_07042026_rsmooth_12/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs12_ngrid2048_thr0p2_15d.fits",
    16: "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_07042026_rsmooth_16/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs16_ngrid2048_thr0p2_15d.fits",
    20: "/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_07042026_rsmooth_20/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000_with_tweb_eigs_rs20_ngrid2048_thr0p2_15d.fits",
}


def radec_z_to_xyz(ra, dec, z):
    zt = np.linspace(0, max(z.max(), 0.6), 4000)
    dt = cosmo.comoving_distance(zt).value
    d = np.interp(z, zt, dt)
    r, dd = np.deg2rad(ra), np.deg2rad(dec)
    return np.vstack([d*np.cos(dd)*np.cos(r), d*np.cos(dd)*np.sin(r), d*np.sin(dd)]).T


def aperture_density(field_pos, query_pos, radius):
    """count of field galaxies within radius of each query galaxy (the sparse sample)."""
    tree = cKDTree(field_pos)
    return np.array(tree.query_ball_point(query_pos, radius, return_length=True))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", nargs=2, type=float, default=[120, 160])
    ap.add_argument("--dec", nargs=2, type=float, default=[14.5, 30.6])
    ap.add_argument("--zr", nargs=2, type=float, default=[0.2, 0.3])
    ap.add_argument("--target-n", type=int, default=120000, help="downsample to BGS-like density")
    ap.add_argument("--apertures-hmpc", type=float, nargs="+", default=[3, 5, 7, 10, 14, 20])
    ap.add_argument("--lambda-th", type=float, default=0.2)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = np.random.default_rng(args.seed); h = cosmo.h

    # footprint from rs7 file
    base = fitsio.read(SCALE_FILES[7], columns=["RA", "DEC", "Z"])
    ra, dec, z = base["RA"], base["DEC"], base["Z"]
    fp = ((ra >= args.ra[0]) & (ra < args.ra[1]) & (dec >= args.dec[0]) & (dec < args.dec[1]) &
          (z >= args.zr[0]) & (z < args.zr[1]))
    fidx = np.where(fp)[0]
    print(f"footprint galaxies (cutsky, dense): {len(fidx)}")
    # downsample to BGS-like density (the sample the model would see)
    if len(fidx) > args.target_n:
        sel = rng.permutation(len(fidx))[: args.target_n]
        fidx = fidx[sel]
    pos = radec_z_to_xyz(ra[fidx], dec[fidx], z[fidx])
    print(f"using {len(fidx)} galaxies; box volume gives mean spacing ~"
          f"{(np.ptp(pos[:,0])*np.ptp(pos[:,1])*np.ptp(pos[:,2])/len(fidx))**(1/3):.1f} Mpc")

    # aperture-density features (sample sees its own sampling -> field=query=downsampled set)
    apsizes = [R/h for R in args.apertures_hmpc]
    dens = np.column_stack([np.log1p(aperture_density(pos, pos, r)) for r in apsizes])
    print(f"aperture-density features at {args.apertures_hmpc} Mpc/h built")

    print("\nscale  std(l1)  frac(l1>th)  | geomR2(l1)  geomR2(l2) | best-single-aperture for l1")
    for s in sorted(SCALE_FILES):
        col = fitsio.read(SCALE_FILES[s], columns=["LAMBDA1", "LAMBDA2"])
        l1 = col["LAMBDA1"][fidx].astype(np.float64)
        l2 = col["LAMBDA2"][fidx].astype(np.float64)
        def r2(X, y):
            m = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08, max_depth=6)
            return r2_score(y, cross_val_predict(m, X, y, cv=args.folds))
        r1 = r2(dens, l1); r2v = r2(dens, l2)
        # which single aperture predicts l1 best?
        singles = [(args.apertures_hmpc[i], r2(dens[:, [i]], l1)) for i in range(dens.shape[1])]
        best = max(singles, key=lambda t: t[1])
        fr = float(np.mean(l1 > args.lambda_th))
        print(f"rs{s:<3d}  {l1.std():.3f}    {fr:.3f}        | {r1:.3f}      {r2v:.3f}     "
              f"| {best[0]:g} Mpc/h (R2={best[1]:.3f})  [" +
              " ".join(f"{a:g}:{rr:.2f}" for a, rr in singles) + "]")

    print("\nReadout: does geomR2(l1) rise with smoothing scale (smoother target = more "
          "learnable from sparse galaxies)? does the best single aperture track the smoothing "
          "scale (scale-matching)? where does predictability plateau = natural label scale.")


if __name__ == "__main__":
    main()
