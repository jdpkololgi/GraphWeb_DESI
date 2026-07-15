#!/usr/bin/env python3
"""Compute DESI BGS ABSMAG_RP1 via the OFFICIAL LSS k+e-correction, for train/inference parity.

Why: the GraphWeb accuracy sprint (Workstream B) wants a luminosity channel. Training uses the
Abacus mock's forward-modelled R_MAG_ABS (no k-correction needed). To DEPLOY that feature on real
DESI we must reproduce an equivalent absolute magnitude — and per the memo we must first verify the
Abacus↔DESI support matches, or the model sees shifted luminosities at inference.

Uses LSS.common_tools.add_ke (Smith+2017 GAMA k-correction via DESI_ke/smith_kcorr, TMR e-correction,
DESI-fiducial distance modulus) -> ABSMAG_RP1 / ABSMAG_RP0. We do NOT reimplement the k-correction.

NOTE: ABSMAG_RP1 is k+e corrected (includes tmr_ecorr). The Abacus R_MAG_ABS convention may differ
(band, rest-frame zref, e-correction) -> the support comparison below is the whole point.

DESI_ke is absent from our LSS checkout, so LSSCODE defaults to a complete NERSC checkout.
Run inside the DESI environment (source /global/common/software/desi/desi_environment.sh main).
"""
from __future__ import annotations
import os, sys, argparse

DEFAULT_LSSCODE = "/global/common/software/desi/users/ioannis"   # has LSS/py/LSS/DESI_ke
WEDGE = "/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits"
FULL = "/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full.dat.fits"
ABACUS = ("/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v0.1/z0.200/"
          "cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lsscode", default=os.environ.get("LSSCODE", DEFAULT_LSSCODE))
    ap.add_argument("--wedge", default=WEDGE)
    ap.add_argument("--full", default=FULL)
    ap.add_argument("--out", default="/pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_wedge_absmag_rp1.fits")
    ap.add_argument("--zlo", type=float, default=0.15)
    ap.add_argument("--zhi", type=float, default=0.55)
    args = ap.parse_args()

    os.environ["LSSCODE"] = args.lsscode
    sys.path.insert(0, f"{args.lsscode}/LSS/py")

    # setuptools>=81 REMOVED the deprecated pkg_resources, but DESI_ke's smith_kcorr.py and
    # findfile.py still `from pkg_resources import resource_filename`. It is a DEAD import there
    # (the only call site, findfile.py:104, is commented out), so a minimal shim satisfies it
    # without patching the shared read-only LSS checkout or polluting the DESI environment.
    import types, importlib.util
    if "pkg_resources" not in sys.modules:
        _shim = types.ModuleType("pkg_resources")

        def _resource_filename(package, resource):        # functional, in case it is ever called
            spec = importlib.util.find_spec(package)
            base = (os.path.dirname(spec.origin) if spec and spec.origin
                    else os.environ.get("CODE_ROOT", "."))
            return os.path.join(base, resource)

        _shim.resource_filename = _resource_filename
        sys.modules["pkg_resources"] = _shim
        print("shimmed pkg_resources (setuptools>=81 removed it; DESI_ke still imports it)")

    import numpy as np
    import fitsio
    from astropy.table import Table
    from LSS.common_tools import add_dered_flux, add_ke
    print(f"LSSCODE={args.lsscode}")

    w = fitsio.read(args.wedge, columns=["TARGETID", "Z", "TARGET_RA", "TARGET_DEC"])
    print(f"wedge rows: {len(w):,}")
    p = fitsio.read(args.full, columns=["TARGETID", "FLUX_G", "FLUX_R",
                                        "MW_TRANSMISSION_G", "MW_TRANSMISSION_R"])
    print(f"LSS full BGS rows: {len(p):,}")

    # join photometry onto the wedge by TARGETID
    order = np.argsort(p["TARGETID"])
    pt = p["TARGETID"][order]
    pos = np.searchsorted(pt, w["TARGETID"])
    pos = np.clip(pos, 0, len(pt) - 1)
    ok = pt[pos] == w["TARGETID"]
    print(f"photometry join: {ok.sum():,}/{len(w):,} = {100*ok.mean():.2f}%")
    src = order[pos[ok]]

    dat = Table()
    dat["TARGETID"] = w["TARGETID"][ok]
    dat["Z"] = w["Z"][ok].astype(float)
    for c in ["FLUX_G", "FLUX_R", "MW_TRANSMISSION_G", "MW_TRANSMISSION_R"]:
        dat[c] = p[c][src].astype(float)

    dat = add_dered_flux(dat, fcols=["G", "R"])          # -> flux_g_dered, flux_r_dered
    good = (dat["flux_g_dered"] > 0) & (dat["flux_r_dered"] > 0) & (dat["Z"] > 0)
    dat = dat[good]
    print(f"after positive-flux cut: {len(dat):,}")
    dat = add_ke(dat, zcol="Z")                          # -> ABSMAG_RP1/RP0, REST_GMR_0P1, KCORR_*
    dat["G_R_OBS"] = (22.5 - 2.5*np.log10(dat["flux_g_dered"])) - (22.5 - 2.5*np.log10(dat["flux_r_dered"]))
    dat.write(args.out, overwrite=True)
    print(f"saved -> {args.out}")

    # ---- SUPPORT COMPARISON: DESI ABSMAG_RP1 vs Abacus R_MAG_ABS (the memo's gate) ----
    z = np.asarray(dat["Z"]); rp1 = np.asarray(dat["ABSMAG_RP1"]); gro_d = np.asarray(dat["G_R_OBS"])
    m = (z >= args.zlo) & (z < args.zhi) & np.isfinite(rp1) & (rp1 < 0) & (rp1 > -30)
    a = fitsio.read(ABACUS, columns=["RA", "DEC", "Z", "R_MAG_ABS", "R_MAG_APP", "G_R_OBS"])
    am = ((a["RA"] >= 118) & (a["RA"] < 162) & (a["DEC"] >= 12.5) & (a["DEC"] < 32.6) &
          (a["Z"] >= args.zlo) & (a["Z"] < args.zhi) & (a["R_MAG_APP"] < 19.5))
    ar = a["R_MAG_ABS"][am].astype(float); ag = a["G_R_OBS"][am].astype(float)
    print(f"\n=== SUPPORT COMPARISON (z {args.zlo}-{args.zhi}) ===")
    print(f"{'quantity':34s} {'n':>9s} {'mean':>8s} {'std':>7s} {'p5':>8s} {'p50':>8s} {'p95':>8s}")
    for nm, v in [("DESI ABSMAG_RP1 (k+e)", rp1[m]), ("Abacus R_MAG_ABS", ar),
                  ("DESI G_R_OBS", gro_d[m]), ("Abacus G_R_OBS", ag)]:
        print(f"{nm:34s} {len(v):9,d} {v.mean():8.3f} {v.std():7.3f} "
              f"{np.percentile(v,5):8.3f} {np.percentile(v,50):8.3f} {np.percentile(v,95):8.3f}")
    print(f"\n  M_abs median OFFSET (Abacus - DESI) = {np.median(ar)-np.median(rp1[m]):+.3f} mag")
    print(f"  M_abs width ratio (Abacus/DESI)      = {ar.std()/rp1[m].std():.3f}")
    print(f"  G_R_OBS median OFFSET (Abacus-DESI)  = {np.median(ag)-np.median(gro_d[m]):+.3f}")
    print("\nGATE: |offset| small and width ratio ~1 => luminosity channel is train/inference safe.")
    print("      Large offset => recalibrate to the DESI convention, or fall back to G_R_OBS only.")


if __name__ == "__main__":
    main()
