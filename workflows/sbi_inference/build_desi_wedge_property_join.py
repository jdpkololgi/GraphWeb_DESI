#!/usr/bin/env python3
"""Join NPE-inferred cosmic-web environment to loa FastSpecFit galaxy properties.

Property->environment closure test, build step ("match galaxy properties to their
TARGETIDs in the wedge"). Produces a tidy per-galaxy table that pairs, for every
galaxy in the DESI LOA wedge:

  * the NPE-inferred T-web environment (hard class, class probabilities, the three
    Hessian/tidal eigenvalue posterior means+stds, and E[trace] = sum of the means),
  * the loa FastSpecFit physical properties (stellar mass, SFR, Dn4000, rest-frame
    SDSS g and r absolute magnitudes), joined by TARGETID.

No inference is re-run: everything is recovered from the existing on-disk products.

Provenance chain (see plan / SCIENCE_LOG 2026-06-19):
  preds.npz["global_node_id"]  -- row index into the source BGS catalog
      |  TARGETID = source_catalog["TARGETID"][global_node_id]
      v
  loa FastSpecFit SPECPHOT (TARGETID -> LOGMSTAR, SFR, DN4000, ABSMAG01_SDSS_G/R)

Caveats baked in (and reported in join_report.json):
  * ~1,252 TARGETIDs appear twice because the full graph is split into Galactic
    N/S hemispheres; each copy gets its own neighbourhood -> its own prediction. We
    AVERAGE the posterior products across the duplicate nodes and keep one row
    (is_dup flag set), then recompute the hard class from the averaged classprob.
  * loa FastSpecFit is the SAME release as the wedge, so TARGETIDs match cleanly,
    but fibre-assignment incompleteness is mildly environment-dependent (denser in
    clusters) and can bias *absolute* fractions slightly; the *trend* is robust.
  * FastSpecFit SFR can be 0 (no detected star formation) -> log sSFR = -inf. Such
    galaxies are counted as quenched (quenched fraction) but excluded from the
    sSFR-median panel (has_ssfr=False).

Run on a compute node (I/O-bound column reads of a few large FITS) with the
cosmic_env python; see the plan's run recipe.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np

# Repo root on path for shared.config_paths (script lives in workflows/sbi_inference/).
_REPO = Path(__file__).resolve().parents[2]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))
from shared.config_paths import (  # noqa: E402
    DESI_FASTSPEC_CATALOGS_DIR,
    GRAPHWEB_SCRATCH_ROOT,
)

CLASS_ORDER = ["void", "wall", "filament", "cluster"]  # ascending lambda1<=lambda2<=lambda3
QUENCHED_LOG_SSFR = -11.0  # log10(sSFR / yr^-1) threshold; below = quenched

# FastSpecFit columns we need (HDU-specific). Mirrors load_catalog.py groupings.
SPECPHOT_COLS = [
    "TARGETID",
    "LOGMSTAR",
    "LOGMSTAR_IVAR",
    "SFR",
    "SFR_IVAR",
    "DN4000",
    "ABSMAG01_SDSS_G",
    "ABSMAG01_SDSS_R",
]
METADATA_COLS = ["TARGETID", "RA", "DEC", "Z"]  # for the join cross-check


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    flx = f"{GRAPHWEB_SCRATCH_ROOT}/flowjax_inference_outputs"
    p.add_argument(
        "--preds",
        type=Path,
        default=Path(f"{flx}/desi_wedge_flowjax_linear_si/desi_wedge_flowjax_preds.npz"),
        help="DESI wedge NPE inference output npz (default: the production SI run).",
    )
    p.add_argument(
        "--source-catalog",
        type=Path,
        default=Path(f"{GRAPHWEB_SCRATCH_ROOT}/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits"),
        help="BGS FITS the graph was built from; its TARGETID column is indexed by global_node_id.",
    )
    p.add_argument(
        "--fastspec-dir",
        type=Path,
        default=Path(DESI_FASTSPEC_CATALOGS_DIR),
        help="loa FastSpecFit catalogs dir (healpix main-bright files).",
    )
    p.add_argument(
        "--fastspec-glob",
        default="fastspec-loa-main-bright-nside1-hp*.fits",
        help="Glob (within --fastspec-dir) for the BGS-bright healpix FastSpecFit files.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output parquet path (default: <preds dir>/desi_wedge_env_props.parquet).",
    )
    return p.parse_args()


def _native(a: np.ndarray) -> np.ndarray:
    """fitsio returns big-endian columns; make them native for pandas/parquet."""
    return np.asarray(a).astype(np.asarray(a).dtype.newbyteorder("="), copy=False)


def main() -> None:
    args = parse_args()
    import fitsio  # local import so --help works without the dep
    import pandas as pd

    preds_path = args.preds.expanduser().resolve()
    out_path = (args.out or (preds_path.parent / "desi_wedge_env_props.parquet")).expanduser().resolve()
    report_path = out_path.with_name("join_report.json")
    print(f"[build] preds            : {preds_path}", flush=True)
    print(f"[build] source catalog   : {args.source_catalog}", flush=True)
    print(f"[build] fastspec dir     : {args.fastspec_dir}", flush=True)

    # --- 1. inference products ------------------------------------------------
    preds = np.load(preds_path)
    gid = preds["global_node_id"].astype(np.int64)
    classprob = np.asarray(preds["classprob"], np.float64)        # [N,4] void..cluster
    lam_mean = np.asarray(preds["lambda_mean"], np.float64)       # [N,3] ascending
    lam_std = np.asarray(preds["lambda_std"], np.float64)
    p_exceed = np.asarray(preds["p_exceed"], np.float64)          # [N,3]
    ra, dec, z = (np.asarray(preds[k], np.float64) for k in ("ra", "dec", "z"))
    n_raw = gid.size
    print(f"[build] inference nodes  : {n_raw:,}", flush=True)

    # --- 2. recover TARGETID via the source-catalog row index -----------------
    all_tid = _native(fitsio.read(str(args.source_catalog), columns=["TARGETID"])["TARGETID"]).astype(np.int64)
    if gid.max() >= all_tid.size:
        raise SystemExit(f"global_node_id max {gid.max()} >= catalog rows {all_tid.size}: wrong source catalog?")
    tid = all_tid[gid]

    df = pd.DataFrame(
        {
            "TARGETID": tid,
            "global_node_id": gid,
            "ra": ra, "dec": dec, "z": z,
            "P_void": classprob[:, 0], "P_wall": classprob[:, 1],
            "P_filament": classprob[:, 2], "P_cluster": classprob[:, 3],
            "lambda_1_mean": lam_mean[:, 0], "lambda_2_mean": lam_mean[:, 1], "lambda_3_mean": lam_mean[:, 2],
            "lambda_1_std": lam_std[:, 0], "lambda_2_std": lam_std[:, 1], "lambda_3_std": lam_std[:, 2],
            "p_exceed_1": p_exceed[:, 0], "p_exceed_2": p_exceed[:, 1], "p_exceed_3": p_exceed[:, 2],
        }
    )

    # --- 3. de-duplicate hemisphere copies: average posterior products --------
    counts = df["TARGETID"].value_counts()
    dup_tids = set(counts[counts > 1].index)
    df["is_dup"] = df["TARGETID"].isin(dup_tids)
    mean_cols = [c for c in df.columns if c.startswith(("P_", "lambda_", "p_exceed_"))]
    agg = {c: "mean" for c in mean_cols}
    agg.update({"ra": "first", "dec": "first", "z": "first",
                "global_node_id": "first", "is_dup": "first"})
    df = df.groupby("TARGETID", as_index=False).agg(agg)
    n_unique = len(df)
    print(f"[build] unique targets   : {n_unique:,} ({len(dup_tids):,} were duplicated x2)", flush=True)

    # Environment summaries derived from the (possibly averaged) posterior.
    P = df[["P_void", "P_wall", "P_filament", "P_cluster"]].to_numpy()
    df["hard_class"] = P.argmax(axis=1).astype(np.int8)
    df["trace_lambda"] = df[["lambda_1_mean", "lambda_2_mean", "lambda_3_mean"]].sum(axis=1)

    # --- 4. join loa FastSpecFit by TARGETID (scheme-agnostic scan) -----------
    fs_files = sorted(glob.glob(str(args.fastspec_dir / args.fastspec_glob)))
    if not fs_files:
        raise SystemExit(f"No FastSpecFit files match {args.fastspec_dir/args.fastspec_glob}")
    wedge_tids = df["TARGETID"].to_numpy(np.int64)
    sp_parts, md_parts = [], []
    for f in fs_files:
        fcol = _native(fitsio.read(f, ext="SPECPHOT", columns=["TARGETID"])["TARGETID"]).astype(np.int64)
        mask = np.isin(fcol, wedge_tids, assume_unique=False)
        nmatch = int(mask.sum())
        print(f"[build]   {os.path.basename(f)}: {nmatch:,} matches", flush=True)
        if nmatch == 0:
            continue
        rows = np.where(mask)[0]
        sp_parts.append(fitsio.read(f, ext="SPECPHOT", rows=rows, columns=SPECPHOT_COLS))
        md_parts.append(fitsio.read(f, ext="METADATA", rows=rows, columns=METADATA_COLS))

    if not sp_parts:
        raise SystemExit("No wedge TARGETIDs found in any FastSpecFit file — wrong release/glob?")
    sp = np.concatenate(sp_parts)
    md = np.concatenate(md_parts)
    fs_df = pd.DataFrame({c: _native(sp[c]).astype(np.float64) if c != "TARGETID"
                          else _native(sp[c]).astype(np.int64) for c in SPECPHOT_COLS})
    fs_df["RA_fs"] = _native(md["RA"]).astype(np.float64)
    fs_df["DEC_fs"] = _native(md["DEC"]).astype(np.float64)
    fs_df["Z_fs"] = _native(md["Z"]).astype(np.float64)
    fs_df = fs_df.drop_duplicates("TARGETID", keep="first")
    out = df.merge(fs_df, on="TARGETID", how="left")

    # --- 5. derived galaxy properties + quality flags -------------------------
    logm = out["LOGMSTAR"].to_numpy()
    sfr = out["SFR"].to_numpy()
    out["gr"] = out["ABSMAG01_SDSS_G"] - out["ABSMAG01_SDSS_R"]
    with np.errstate(divide="ignore", invalid="ignore"):
        out["log_sSFR"] = np.log10(sfr) - logm  # log10(SFR/M*) in yr^-1
    out["has_props"] = np.isfinite(logm) & (out["LOGMSTAR_IVAR"].to_numpy() > 0)
    out["has_ssfr"] = out["has_props"].to_numpy() & (sfr > 0) & np.isfinite(out["log_sSFR"].to_numpy())
    # SFR<=0 -> genuinely quenched (sSFR = -inf); counted in quenched fraction.
    quenched = out["has_props"].to_numpy() & ((sfr <= 0) | (out["log_sSFR"].to_numpy() < QUENCHED_LOG_SSFR))
    out["quenched"] = quenched

    # --- 6. join cross-check (preds RA/Dec vs FastSpecFit RA/Dec) -------------
    m = out["has_props"].to_numpy()
    if m.any():
        dra = (out["ra"].to_numpy()[m] - out["RA_fs"].to_numpy()[m]) * np.cos(np.radians(out["dec"].to_numpy()[m]))
        ddec = out["dec"].to_numpy()[m] - out["DEC_fs"].to_numpy()[m]
        sep_arcsec = np.sqrt(dra ** 2 + ddec ** 2) * 3600.0
        max_sep = float(np.nanmax(sep_arcsec))
    else:
        max_sep = float("nan")

    # --- 7. write outputs -----------------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        out.to_parquet(out_path, index=False)
        written = out_path
    except Exception as e:  # pyarrow/fastparquet missing -> csv fallback
        written = out_path.with_suffix(".csv")
        out.to_csv(written, index=False)
        print(f"[build] parquet unavailable ({e}); wrote CSV instead.", flush=True)
    print(f"[build] wrote table      : {written}  ({len(out):,} rows)", flush=True)

    matched = int(out["has_props"].sum())
    report = {
        "preds": str(preds_path),
        "source_catalog": str(args.source_catalog),
        "fastspec_dir": str(args.fastspec_dir),
        "fastspec_files": [os.path.basename(f) for f in fs_files],
        "quenched_log_ssfr_threshold": QUENCHED_LOG_SSFR,
        "n_inference_nodes": n_raw,
        "n_unique_targets": n_unique,
        "n_duplicate_targets": len(dup_tids),
        "n_matched_fastspec": matched,
        "frac_matched_fastspec": matched / n_unique if n_unique else 0.0,
        "n_valid_ssfr": int(out["has_ssfr"].sum()),
        "frac_valid_ssfr": float(out["has_ssfr"].mean()),
        "radec_crosscheck_max_arcsec": max_sep,
        "class_counts_all": {CLASS_ORDER[k]: int((out["hard_class"] == k).sum()) for k in range(4)},
        "class_counts_matched": {
            CLASS_ORDER[k]: int(((out["hard_class"] == k) & out["has_props"]).sum()) for k in range(4)
        },
        "quenched_fraction_overall": float(out.loc[out["has_props"], "quenched"].mean()) if matched else None,
        "table_path": str(written),
    }
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"[build] wrote report     : {report_path}", flush=True)
    print(f"[build] FastSpecFit match : {matched:,}/{n_unique:,} "
          f"({100*report['frac_matched_fastspec']:.1f}%); valid sSFR "
          f"{report['n_valid_ssfr']:,} ({100*report['frac_valid_ssfr']:.1f}%)", flush=True)
    print(f"[build] RA/Dec cross-check max sep: {max_sep:.3f} arcsec "
          f"({'OK' if max_sep < 1.0 else 'WARN — check join!'})", flush=True)


if __name__ == "__main__":
    main()
