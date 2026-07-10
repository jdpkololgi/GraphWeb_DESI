#!/usr/bin/env python3
"""Approach A: re-join CIGALE (HZ) stellar mass + SFR onto the EXISTING wedge
posteriors by TARGETID, replacing the problematic FastSpecFit SFRs.

The environment inference (positions, graph, posteriors) is unchanged and stays
parity-matched to the mock (DELTACHI2>=25 harmonization) — we only swap the
galaxy PROPERTIES used in the science plots. Writes a parquet with the SAME
schema as desi_wedge_env_props.parquet so the existing plot scripts work by
repointing, but LOGMSTAR/SFR/log_sSFR/quenched now come from CIGALE-HZ CG_15.

Coverage: ~90% of the 111,503 wedge galaxies match the goodPhoto CIGALE catalogue.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import fitsio

EXIST = Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
             "desi_wedge_flowjax_linear_si/desi_wedge_env_props.parquet")
CIG = Path("/global/cfs/cdirs/desi/users/manasvee/1_prepped_data/DESI/loa/"
           "desi_loa_fsf_bgs_scnd_quality_ZallPix_cigaleHu_primaryZ_goodPhoto.fits")
OUT_DIR = Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/"
               "desi_wedge_cigale_hz")
VAR = "15"  # CIGALE config: CG_15 (fiducial; ~identical to CG_5)
QUENCH_LOGSSFR = -11.0


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.read_parquet(EXIST)
    print(f"existing wedge rows: {len(df):,}")

    # CIGALE side: wedge cut, dedup by TARGETID (keep first)
    c = fitsio.read(CIG, columns=["TARGETID", "RA", "DEC", "Z_fsf",
                                  f"SFR_CG_{VAR}", f"MASS_CG_{VAR}",
                                  "SFR_CG_5", "MASS_CG_5"])
    w = ((c["RA"] >= 120) & (c["RA"] < 160) & (c["DEC"] >= 14.5)
         & (c["DEC"] < 30.6) & (c["Z_fsf"] >= 0.2) & (c["Z_fsf"] < 0.3))
    c = c[w]
    nat = lambda a: np.asarray(a).astype(a.dtype.newbyteorder("="))  # FITS big-endian -> native
    ct = pd.DataFrame({
        "TARGETID": nat(c["TARGETID"]).astype(np.int64),
        "SFR_CG": nat(c[f"SFR_CG_{VAR}"]).astype(np.float64),
        "MASS_CG": nat(c[f"MASS_CG_{VAR}"]).astype(np.float64),
        "SFR_CG5": nat(c["SFR_CG_5"]).astype(np.float64),
        "MASS_CG5": nat(c["MASS_CG_5"]).astype(np.float64),
    }).drop_duplicates("TARGETID", keep="first")
    print(f"CIGALE wedge unique TARGETIDs: {len(ct):,}")

    m = df.merge(ct, on="TARGETID", how="left")
    matched = m["MASS_CG"].notna().to_numpy()
    print(f"matched: {matched.sum():,} / {len(m):,}  ({100*matched.mean():.1f}%)")

    # CIGALE-derived properties overwrite the FastSpecFit ones (same column names
    # so the plot scripts run unchanged); provenance kept in *_CG columns.
    mass = m["MASS_CG"].to_numpy(); sfr = m["SFR_CG"].to_numpy()
    good = matched & (mass > 0) & (sfr > 0) & np.isfinite(mass) & np.isfinite(sfr)
    m["LOGMSTAR"] = np.where(good, np.log10(np.where(mass > 0, mass, np.nan)), np.nan)
    m["SFR"] = np.where(good, sfr, np.nan)
    m["log_sSFR"] = np.where(good, np.log10(np.where(sfr > 0, sfr, np.nan))
                             - m["LOGMSTAR"], np.nan)
    m["quenched"] = m["log_sSFR"] < QUENCH_LOGSSFR
    m["has_ssfr"] = good
    m["has_props"] = matched
    # keep gr / DN4000 / ABSMAG from the existing photometry (unchanged)
    m["sfr_source"] = f"CIGALE_HZ_CG{VAR}"

    out = OUT_DIR / "desi_wedge_env_props.parquet"
    m.to_parquet(out, index=False)
    print(f"WROTE {out}  ({len(m):,} rows, {good.sum():,} with CIGALE sSFR)")
    # sanity: sSFR bimodality
    s = m.loc[good, "log_sSFR"].to_numpy(); s = s[(s > -14) & (s < -8)]
    h, e = np.histogram(s, bins=np.arange(-14, -8, 0.15)); cc = 0.5 * (e[:-1] + e[1:])
    pk = [round(cc[i], 2) for i in range(1, len(h) - 1)
          if h[i] > h[i-1] and h[i] >= h[i+1] and h[i] > h.max() * 0.12]
    print("CIGALE log sSFR peaks:", pk, " median logM*:",
          round(float(np.nanmedian(m.loc[good, 'LOGMSTAR'])), 2))


if __name__ == "__main__":
    main()
