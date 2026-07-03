#!/usr/bin/env python3
"""A2 — clean n(z) comparison, path1 mock parent vs real DESI BGS (roadmap v2, Track 2).

Same wedge RA/Dec box for both samples => same solid angle => shell counts compare
directly; also reports comoving number density n(z) [Planck18]. Mock phantoms
(sentinel z~0.5898, the injection bug) are counted and excluded EXPLICITLY so the
curve is the clean observed-only n(z) without waiting for the A1 regeneration.

Outputs: per-shell table + themed figure + JSON to a new scratch dir. CPU only.
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from pathlib import Path
import numpy as np
import fitsio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18 as cosmo

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

MOCK = ("/pscratch/sd/d/dkololgi/abacus/mocks_with_eigs_05062026_rsmooth_7/"
        "mock_bgs_maglim_path1_fiberassign_graph_ready_with_tweb_eigs_rs7_ngrid2048_thr0p2_halo_xcom.fits")
DESI = "/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits"
SENTINEL = (0.585, 0.595)


def shell_counts(z, edges):
    return np.histogram(z, bins=edges)[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ra", nargs=2, type=float, default=[120, 160])
    ap.add_argument("--dec", nargs=2, type=float, default=[14.5, 30.6])
    ap.add_argument("--zmin", type=float, default=0.05)
    ap.add_argument("--zmax", type=float, default=0.55)
    ap.add_argument("--dz", type=float, default=0.01)
    ap.add_argument("--out-dir", type=Path,
                    default=Path("/pscratch/sd/d/dkololgi/abacus/nz_comparison_20260703"))
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    apply_style()

    def load(path, zcol="Z"):
        names = [c.upper() for c in fitsio.FITS(path)[1].get_colnames()]
        ra_c = "RA" if "RA" in names else "TARGET_RA"
        dec_c = "DEC" if "DEC" in names else "TARGET_DEC"
        t = fitsio.read(path, columns=[ra_c, dec_c, zcol])
        sel = ((t[ra_c] >= args.ra[0]) & (t[ra_c] < args.ra[1]) &
               (t[dec_c] >= args.dec[0]) & (t[dec_c] < args.dec[1]))
        return np.asarray(t[zcol][sel], np.float64)

    zm = load(MOCK)
    zd = load(DESI)
    phantom = (zm >= SENTINEL[0]) & (zm < SENTINEL[1])
    print(f"mock wedge-box rows: {len(zm)}  (phantoms in sentinel window {SENTINEL}: "
          f"{phantom.sum()} = {phantom.mean():.1%}, EXCLUDED)")
    zm = zm[~phantom]
    print(f"DESI wedge-box rows: {len(zd)}")

    edges = np.arange(args.zmin, args.zmax + args.dz / 2, args.dz)
    mid = 0.5 * (edges[:-1] + edges[1:])
    cm = shell_counts(zm, edges).astype(float)
    cd = shell_counts(zd, edges).astype(float)

    # solid angle of the RA/Dec box (same for both samples)
    dra = np.deg2rad(args.ra[1] - args.ra[0])
    sind = np.sin(np.deg2rad(args.dec[1])) - np.sin(np.deg2rad(args.dec[0]))
    omega = dra * sind                                   # steradians
    dcom = cosmo.comoving_distance(edges).value          # Mpc
    vshell = (omega / 3.0) * (dcom[1:] ** 3 - dcom[:-1] ** 3)
    nm = cm / vshell; nd = cd / vshell
    ratio = np.divide(cm, cd, out=np.full_like(cm, np.nan), where=cd > 0)

    z23 = (mid >= 0.2) & (mid < 0.3)
    print(f"\nwedge z0.2-0.3 totals: mock {cm[z23].sum():.0f} vs DESI {cd[z23].sum():.0f} "
          f"(mock/DESI = {cm[z23].sum()/cd[z23].sum():.3f})")
    print("z_mid   N_mock   N_DESI   mock/DESI   n_mock[Mpc^-3]   n_DESI")
    for i in range(len(mid)):
        if cd[i] > 0 or cm[i] > 0:
            print(f"{mid[i]:.3f}  {cm[i]:7.0f}  {cd[i]:7.0f}   {ratio[i]:8.3f}   "
                  f"{nm[i]:.3e}     {nd[i]:.3e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.plot(mid, nd, "-", lw=2, color=ACCENT_COLORS["magenta"], label="DESI BGS")
    ax1.plot(mid, nm, "-", lw=2, color=ACCENT_COLORS["blue"], label="path1 mock (clean)")
    ax1.axvspan(0.2, 0.3, alpha=0.12, color="#9a9a93", label="current wedge")
    ax1.set_xlabel("z"); ax1.set_ylabel(r"$n(z)$ [Mpc$^{-3}$]"); ax1.set_yscale("log")
    ax1.legend(); ax1.set_title("comoving number density (same footprint box)")
    ax2.plot(mid, ratio, "-", lw=2, color=ACCENT_COLORS["red"])
    ax2.axhline(1.0, ls=":", color="#9a9a93"); ax2.axvspan(0.2, 0.3, alpha=0.12, color="#9a9a93")
    ax2.set_xlabel("z"); ax2.set_ylabel("mock / DESI"); ax2.set_ylim(0, 2)
    ax2.set_title("shell-count ratio")
    fig.suptitle("A2: n(z), path1 fiberassign mock vs DESI BGS bright (wedge box)")
    fig.savefig(args.out_dir / "nz_mock_vs_desi.png", bbox_inches="tight", dpi=200)

    json.dump({"z_mid": mid.tolist(), "n_mock": nm.tolist(), "n_desi": nd.tolist(),
               "count_mock": cm.tolist(), "count_desi": cd.tolist(),
               "ratio": ratio.tolist(), "omega_sr": omega,
               "phantoms_excluded": int(phantom.sum())},
              open(args.out_dir / "nz_mock_vs_desi.json", "w"), indent=2)
    print(f"\nSaved: {args.out_dir}/nz_mock_vs_desi.png + .json")


if __name__ == "__main__":
    main()
