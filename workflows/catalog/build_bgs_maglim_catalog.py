#!/usr/bin/env python3
"""Build a magnitude-limited DESI BGS catalog from the LOA zall file.

This produces a single FITS table filtered by:
- ZWARN == 0
- DELTACHI2 >= 25
- SPECTYPE == GALAXY
- BGS_TARGET has at least one BGS_BRIGHT bit set (north/south/unsplit), unless --no-bright-only

No redshift cuts, no stellar-mass cuts.

Default input is the LOA zcat:
  /global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits

Output is written as a FITS binary table using streaming chunks to avoid
loading the full zall catalog into memory.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import fitsio
import numpy as np
from desitarget.targetmask import bgs_mask

BRIGHT_BITS = (
    bgs_mask.BGS_BRIGHT
    | bgs_mask.BGS_BRIGHT_NORTH
    | bgs_mask.BGS_BRIGHT_SOUTH
)

DEFAULT_ZALL = "/global/cfs/cdirs/desi/public/dr2/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits"

# Keep this minimal; you can always join extra VAC columns later.
DEFAULT_COLS = (
    "TARGETID",
    "SURVEY",
    "PROGRAM",
    "TARGET_RA",
    "TARGET_DEC",
    "Z",
    "ZWARN",
    "DELTACHI2",
    "SPECTYPE",
    "BGS_TARGET",
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--zall-path", type=str, default=DEFAULT_ZALL)
    p.add_argument("--out-path", type=Path, required=True)
    p.add_argument("--chunk-rows", type=int, default=2_000_000)
    p.add_argument("--overwrite", action="store_true")
    p.add_argument(
        "--bright-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Require BGS_BRIGHT (or north/south bright) in BGS_TARGET (default: True).",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    zall_path = Path(args.zall_path).expanduser().resolve()
    if not zall_path.exists():
        raise FileNotFoundError(f"zall not found: {zall_path}")

    out = args.out_path.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        if args.overwrite:
            out.unlink()
        else:
            raise FileExistsError(f"Output exists: {out} (use --overwrite)")

    print(f"Input zall: {zall_path}")
    print(f"Output: {out}")
    print(f"bright_only: {args.bright_only}")

    f = fitsio.FITS(str(zall_path), "r")
    try:
        hdu = f["ZCATALOG"] if "ZCATALOG" in f else f[1]
        nrows = int(hdu.get_nrows())
        print(f"Input rows: {nrows:,}")

        fout = fitsio.FITS(str(out), "rw", clobber=True)
        try:
            wrote = 0
            kept = 0
            first = True
            for start in range(0, nrows, int(args.chunk_rows)):
                stop = min(start + int(args.chunk_rows), nrows)
                tab = hdu.read(rows=range(start, stop), columns=list(DEFAULT_COLS))

                m = (tab["ZWARN"] == 0) & (tab["DELTACHI2"] >= 25.0) & (tab["SPECTYPE"] == "GALAXY")
                if args.bright_only:
                    m &= (tab["BGS_TARGET"] & BRIGHT_BITS) != 0
                else:
                    m &= tab["BGS_TARGET"] != 0

                sub = tab[m]
                wrote += stop - start
                kept += int(sub.size)

                if sub.size:
                    if first:
                        fout.write(sub)
                        first = False
                    else:
                        fout[1].append(sub)

                if start == 0 or (start // int(args.chunk_rows) + 1) % 5 == 0 or stop == nrows:
                    frac = kept / max(wrote, 1)
                    print(f"  scanned {stop:,}/{nrows:,} rows, kept={kept:,} (frac={frac:.4f})")

        finally:
            fout.close()

    finally:
        f.close()

    with fitsio.FITS(str(out), "r") as fcheck:
        n_out = int(fcheck[1].get_nrows())
    print(f"Done. Output rows: {n_out:,}")


if __name__ == "__main__":
    main()
