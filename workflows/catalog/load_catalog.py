import argparse
import os
from glob import glob
import sys
from pathlib import Path

import fitsio
import numpy as np
from astropy.table import Table, hstack, join, vstack

# Fallback DESI Python paths for environments where DESI packages are not on
# sys.path by default (e.g. custom conda env without full desienv activation).
DESI_FALLBACK_SITE_PACKAGES = (
    "/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/"
    "code/desiutil/3.5.0/lib/python3.10/site-packages",
    "/global/common/software/desi/perlmutter/desiconda/20240425-2.2.0/"
    "code/desimodel/0.19.3/lib/python3.10/site-packages",
)


def _import_radec2pix():
    try:
        from desimodel.footprint import radec2pix as _radec2pix

        return _radec2pix
    except ImportError:
        for site_path in DESI_FALLBACK_SITE_PACKAGES:
            if os.path.isdir(site_path) and site_path not in sys.path:
                sys.path.insert(0, site_path)
        try:
            from desimodel.footprint import radec2pix as _radec2pix

            return _radec2pix
        except ImportError:
            return None


radec2pix = _import_radec2pix()

# Allow canonical workflow scripts to resolve repo-root modules after reorganization.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config_paths import DESI_FASTSPEC_CATALOGS_DIR, DESI_TRACTORPHOT_DIR, DESI_ZCAT_FILE

# Workflow status: ACTIVE (low-z catalog assembly)


FASTSPEC_CATALOGS = [
    "fastspec-loa-cmx-other.fits",
    "fastspec-loa-main-backup.fits",
    "fastspec-loa-main-bright-nside1-hp00.fits",
    "fastspec-loa-main-bright-nside1-hp01.fits",
    "fastspec-loa-main-bright-nside1-hp02.fits",
    "fastspec-loa-main-bright-nside1-hp03.fits",
    "fastspec-loa-main-bright-nside1-hp04.fits",
    "fastspec-loa-main-bright-nside1-hp05.fits",
    "fastspec-loa-main-bright-nside1-hp06.fits",
    "fastspec-loa-main-bright-nside1-hp07.fits",
    "fastspec-loa-main-bright-nside1-hp08.fits",
    "fastspec-loa-main-bright-nside1-hp09.fits",
    "fastspec-loa-main-bright-nside1-hp10.fits",
    "fastspec-loa-main-bright-nside1-hp11.fits",
    "fastspec-loa-main-dark-nside1-hp00.fits",
    "fastspec-loa-main-dark-nside1-hp01.fits",
    "fastspec-loa-main-dark-nside1-hp02.fits",
    "fastspec-loa-main-dark-nside1-hp03.fits",
    "fastspec-loa-main-dark-nside1-hp04.fits",
    "fastspec-loa-main-dark-nside1-hp05.fits",
    "fastspec-loa-main-dark-nside1-hp06.fits",
    "fastspec-loa-main-dark-nside1-hp07.fits",
    "fastspec-loa-main-dark-nside1-hp08.fits",
    "fastspec-loa-main-dark-nside1-hp09.fits",
    "fastspec-loa-main-dark-nside1-hp10.fits",
    "fastspec-loa-main-dark-nside1-hp11.fits",
    "fastspec-loa-special-backup.fits",
    "fastspec-loa-special-bright.fits",
    "fastspec-loa-special-dark.fits",
    "fastspec-loa-sv1-backup.fits",
    "fastspec-loa-sv1-bright.fits",
    "fastspec-loa-sv1-dark.fits",
    "fastspec-loa-sv1-other.fits",
    "fastspec-loa-sv2-backup.fits",
    "fastspec-loa-sv2-bright.fits",
    "fastspec-loa-sv2-dark.fits",
    "fastspec-loa-sv3-backup.fits",
    "fastspec-loa-sv3-bright.fits",
    "fastspec-loa-sv3-dark.fits",
]

COLS_SELECTION_METADATA = [
    "TARGETID",
    "SURVEY",
    "PROGRAM",
    "RA",
    "DEC",
    "BGS_TARGET",
    "Z",
    "ZWARN",
    "DELTACHI2",
    "SPECTYPE",
    "FLUX_G",
    "FLUX_R",
    "FLUX_Z",
    "FLUX_IVAR_G",
    "FLUX_IVAR_R",
    "FLUX_IVAR_Z",
]

COLS_SELECTION_SPECPHOT = [
    "LOGMSTAR",
    "LOGMSTAR_IVAR",
    "ABSMAG01_SDSS_U",
    "ABSMAG01_SDSS_G",
    "ABSMAG01_SDSS_R",
    "ABSMAG01_SDSS_I",
    "ABSMAG01_SDSS_Z",
    "ABSMAG01_IVAR_SDSS_U",
    "ABSMAG01_IVAR_SDSS_G",
    "ABSMAG01_IVAR_SDSS_R",
    "ABSMAG01_IVAR_SDSS_I",
    "ABSMAG01_IVAR_SDSS_Z",
]

COLS_SELECTION_FASTSPEC = [
    "APERCORR",
    "APERCORR_R",
    "INIT_BALMER_BROAD",
    "INIT_SIGMA_NARROW",
    "INIT_SIGMA_BALMER",
    "HALPHA_AMP",
    "HALPHA_AMP_IVAR",
    "HALPHA_FLUX",
    "HALPHA_FLUX_IVAR",
    "HALPHA_EW",
    "HALPHA_SIGMA",
    "HALPHA_SIGMA_IVAR",
    "HALPHA_BROAD_AMP",
    "HALPHA_BROAD_AMP_IVAR",
    "HALPHA_BROAD_FLUX",
    "HALPHA_BROAD_FLUX_IVAR",
    "HALPHA_BROAD_SIGMA",
    "HALPHA_BROAD_SIGMA_IVAR",
    "HBETA_AMP",
    "HBETA_AMP_IVAR",
    "HBETA_FLUX",
    "HBETA_FLUX_IVAR",
    "HBETA_EW",
    "HBETA_SIGMA",
    "HBETA_SIGMA_IVAR",
    "HBETA_BROAD_AMP",
    "HBETA_BROAD_AMP_IVAR",
    "HBETA_BROAD_FLUX",
    "HBETA_BROAD_FLUX_IVAR",
    "HBETA_BROAD_SIGMA",
    "HBETA_BROAD_SIGMA_IVAR",
    "OIII_5007_AMP",
    "OIII_5007_AMP_IVAR",
    "OIII_5007_FLUX",
    "OIII_5007_FLUX_IVAR",
    "OIII_5007_SIGMA",
    "OIII_5007_SIGMA_IVAR",
    "OII_3726_AMP",
    "OII_3726_AMP_IVAR",
    "OII_3726_FLUX",
    "OII_3726_FLUX_IVAR",
    "OII_3726_SIGMA",
    "OII_3726_SIGMA_IVAR",
    "OII_3729_AMP",
    "OII_3729_AMP_IVAR",
    "OII_3729_FLUX",
    "OII_3729_FLUX_IVAR",
    "OII_3729_SIGMA",
    "OII_3729_SIGMA_IVAR",
    "NII_6548_AMP",
    "NII_6548_AMP_IVAR",
    "NII_6548_FLUX",
    "NII_6548_FLUX_IVAR",
    "NII_6548_SIGMA",
    "NII_6548_SIGMA_IVAR",
    "NII_6584_AMP",
    "NII_6584_AMP_IVAR",
    "NII_6584_FLUX",
    "NII_6584_FLUX_IVAR",
    "NII_6584_SIGMA",
    "NII_6584_SIGMA_IVAR",
    "OIII_5007_CONT",
    "OIII_5007_CONT_IVAR",
    "HALPHA_CONT",
    "HALPHA_CONT_IVAR",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Assemble low-z LOA catalogs into merged FITS products.")
    parser.add_argument("--fastspec-path", default=DESI_FASTSPEC_CATALOGS_DIR)
    parser.add_argument("--zcat-file", default=DESI_ZCAT_FILE)
    parser.add_argument("--vacdir", default=DESI_TRACTORPHOT_DIR)
    parser.add_argument("--specprod", default="loa")
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--z-max", type=float, default=0.06)
    # Catalogs live on pscratch, not in the repo/home (home hit its 40 GiB quota 2026-07-15).
    parser.add_argument("--out-dir", default=os.environ.get(
        "GRAPHWEB_CATALOG_DIR", "/pscratch/sd/d/dkololgi/graphweb_desi/catalogs"))
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args()


def get_selected_rows_metadata(cat, z_min: float, z_max: float) -> np.ndarray:
    return (cat["SPECTYPE"] == "GALAXY") & (cat["Z"] >= z_min) & (cat["Z"] <= z_max)


def read_tractorphot(cat, vacdir: str, specprod: str = "loa", verbose: bool = False) -> Table:
    if radec2pix is None:
        raise ImportError(
            "load_catalog.py requires the DESI package `desimodel` (missing `desimodel.footprint`). "
            "Activate a DESI environment or install desimodel before running catalog assembly."
        )
    tractorphotfiles = glob(os.path.join(vacdir, "tractorphot", f"tractorphot-nside4-hp???-{specprod}.fits"))
    if not tractorphotfiles:
        return Table()

    hdr = fitsio.read_header(tractorphotfiles[0], "TRACTORPHOT")
    tractorphot_nside = hdr["FILENSID"]
    pixels = radec2pix(tractorphot_nside, cat["RA"], cat["DEC"])
    phot = []
    for pixel in sorted(set(pixels)):
        mask = pixel == pixels
        photfile = os.path.join(vacdir, "tractorphot", f"tractorphot-nside4-hp{pixel:03d}-{specprod}.fits")
        if not os.path.isfile(photfile):
            continue
        targetids = fitsio.read(photfile, columns="TARGETID")
        rows = np.where(np.isin(targetids, cat["TARGETID"][mask]))[0]
        if len(rows) == 0:
            continue
        if verbose:
            print(f"Gathering Tractor photometry for {len(rows):,d} object(s) from {photfile}")
        phot.append(Table(fitsio.read(photfile, rows=rows)))

    return vstack(phot) if phot else Table()


def build_selected_rows(fastspec_path: str, z_min: float, z_max: float) -> list[np.ndarray]:
    print("Looking for galaxies matching criteria in fastspecfit files")
    all_rows = []
    matched_counter = 0
    for catalog_name in FASTSPEC_CATALOGS:
        print(f"Processing {catalog_name}")
        cat_meta = fitsio.read(
            os.path.join(fastspec_path, catalog_name),
            "METADATA",
            columns=COLS_SELECTION_METADATA,
        )
        rows_bool = get_selected_rows_metadata(cat_meta, z_min=z_min, z_max=z_max)
        matched_counter += int(np.sum(rows_bool))
        all_rows.append(np.arange(len(rows_bool))[rows_bool])
    print(f"Found {matched_counter:,} galaxies matching criteria")
    print(f"{len(all_rows)} sub-catalogs processed")
    return all_rows


def main() -> None:
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    all_rows = build_selected_rows(args.fastspec_path, args.z_min, args.z_max)

    specphot_tables = []
    fast_tables = []
    meta_tables = []
    for i, catalog_name in enumerate(FASTSPEC_CATALOGS):
        catalog_path = os.path.join(args.fastspec_path, catalog_name)
        specphot_tables.append(Table(fitsio.read(catalog_path, "SPECPHOT", rows=all_rows[i])))
        meta_tables.append(Table(fitsio.read(catalog_path, "METADATA", rows=all_rows[i])))
        fast_tables.append(Table(fitsio.read(catalog_path, "FASTSPEC", rows=all_rows[i])))

    fast_combined = vstack(fast_tables)
    specphot_combined = vstack(specphot_tables)
    meta_combined = vstack(meta_tables)
    meta_combined.keep_columns(COLS_SELECTION_METADATA)
    specphot_combined.keep_columns(COLS_SELECTION_SPECPHOT)
    fast_combined.keep_columns(COLS_SELECTION_FASTSPEC)

    merged_cat = hstack([hstack([meta_combined, specphot_combined]), fast_combined])
    base_path = os.path.join(args.out_dir, "loa-combined-lowz.fits")
    merged_cat.write(base_path, overwrite=True)
    print(f"Finished base extraction: {base_path}")

    print("Matching with redshift catalog")
    zcat = fitsio.read(args.zcat_file, columns=["Z", "SPECTYPE"])
    rows = np.arange(len(zcat))[get_selected_rows_metadata(zcat, z_min=args.z_min, z_max=args.z_max)]
    ztab = Table(fitsio.read(args.zcat_file, rows=rows))
    ztab.keep_columns(
        ["TARGETID", "SURVEY", "PROGRAM", "DESINAME", "MORPHTYPE", "MAIN_PRIMARY", "SV_PRIMARY", "ZCAT_NSPEC", "ZCAT_PRIMARY"]
    )
    cat_wflags = join(merged_cat, ztab, join_type="left", keys=["TARGETID", "SURVEY", "PROGRAM"])
    cat_wflags = cat_wflags[cat_wflags["ZCAT_PRIMARY"] == True]
    zflags_path = os.path.join(args.out_dir, "loa-combined-lowz-zflags.fits")
    cat_wflags.write(zflags_path, overwrite=True)
    print(f"Wrote zflag catalog: {zflags_path}")

    print("Matching with Legacy Survey photometry")
    tractor = read_tractorphot(cat_wflags, vacdir=args.vacdir, specprod=args.specprod, verbose=args.verbose)
    if len(tractor) > 0:
        tractor.keep_columns(
            [
                "TARGETID",
                "TYPE",
                "FRACFLUX_G",
                "FRACFLUX_R",
                "FRACFLUX_Z",
                "SHAPE_R",
                "SHAPE_E1",
                "SHAPE_E2",
                "SHAPE_R_IVAR",
                "SHAPE_E1_IVAR",
                "SHAPE_E2_IVAR",
                "SERSIC",
                "SERSIC_IVAR",
            ]
        )
        cat_withphot = join(cat_wflags, tractor, join_type="left", keys=["TARGETID"])
    else:
        cat_withphot = cat_wflags

    final_path = os.path.join(args.out_dir, "loa-combined-lowz-fastspec-phot.fits")
    cat_withphot.write(final_path, overwrite=True)
    print(f"{len(cat_withphot):,} objects in final catalog")
    print(f"Wrote final catalog: {final_path}")


if __name__ == "__main__":
    main()
