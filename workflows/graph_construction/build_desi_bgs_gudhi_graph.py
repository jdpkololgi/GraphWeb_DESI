#!/usr/bin/env python3
"""Build DESI BGS graph artifacts using Gudhi alpha-complex machinery (hemisphere split).

This intentionally mimics the Abacus `build_abacus_graph.py` behavior:

- Convert (TARGET_RA, TARGET_DEC, Z) to comoving Cartesian (Mpc/h) using Planck18.
- Split Galactic hemispheres (b>0 vs b<=0) and build **one** AlphaComplex per hemisphere.
- Extract undirected edges and merge them into a single edge list using global row indices.

Important:
- This produces **no edges** connecting north<->south, matching the Abacus workflow.
- For a magnitude-limited DESI BGS catalog with ~O(10^7) galaxies, a full Gudhi
  AlphaComplex per hemisphere may be prohibitively expensive in memory/time.
  Use --max-points-per-hemi for smoke tests and incremental development.

Outputs (in --out-dir):
- points_xyz_mpc_h.npy                  : full xyz memmap (float32) for reproducibility
- hemisphere_flag.npy                   : int8 array (1=north, 0=south) aligned to catalog row order
- edges_combined_idx.npy                : undirected edges (u,v) in global row indices (int32, sorted) [legacy]
- tetrahedra_idx.npy                    : tetrahedra vertex indices (T,4) in global row indices (int32) [legacy]
- tetrahedra_volumes.npy                : tetrahedra volumes (T,) in (Mpc/h)^3 (float64) [legacy]
- graph_metadata.json                   : summary and provenance [legacy]

Additionally, Abacus-style outputs (to keep workflows identical):
- <prefix>_edges_combined_idx.npy       : same as Abacus (int32, sorted)
- <prefix>_tetrahedra_idx.npy           : same as Abacus (int32)
- <prefix>_tetrahedra_volumes.npy       : same as Abacus (float64)
- <prefix>_points_xyz.npy               : xyz only (float64; Abacus uses float64)
- <prefix>_metadata.json                : same metadata layout keys as Abacus
"""

from __future__ import annotations

import argparse
import json
import os
import time
import threading
from pathlib import Path

import numpy as np
import fitsio
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18 as cosmo
import astropy.units as u


def _fmt_int(x: int) -> str:
    return f"{int(x):,}"


def _fmt_sec(dt: float) -> str:
    if dt < 120:
        return f"{dt:.1f}s"
    if dt < 7200:
        return f"{dt/60.0:.1f}m"
    return f"{dt/3600.0:.2f}h"


def _print_kv(k: str, v: object) -> None:
    print(f"{k:<24} {v}", flush=True)


def _extract_edges(simplex_tree, alpha_sq: float, *, progress_every: int = 0) -> set[tuple[int, int]]:
    edges: set[tuple[int, int]] = set()
    count = 0
    for simplex, filtration in simplex_tree.get_filtration():
        count += 1
        if progress_every > 0 and count % progress_every == 0:
            print(f"    processed {_fmt_int(count)} simplices...", flush=True)
        if filtration > alpha_sq:
            continue
        if len(simplex) == 2:
            u, v = int(simplex[0]), int(simplex[1])
            if u > v:
                u, v = v, u
            edges.add((u, v))
    return edges


def _tetra_volume(coords: np.ndarray) -> float:
    """Return tetrahedron volume from a (4, 3) coordinate array."""
    a = coords[1] - coords[0]
    b = coords[2] - coords[0]
    c = coords[3] - coords[0]
    return float(abs(np.dot(a, np.cross(b, c))) / 6.0)


def _extract_edges_and_tets(simplex_tree, pts_xyz: np.ndarray, alpha_sq: float, *, progress_every: int = 0):
    """Extract filtered edges + tetrahedra (and volumes) from a simplex tree."""
    edges: set[tuple[int, int]] = set()
    tetrahedra: list[list[int]] = []
    volumes: list[float] = []
    count = 0
    for simplex, filtration in simplex_tree.get_filtration():
        count += 1
        if progress_every > 0 and count % progress_every == 0:
            print(f"    processed {_fmt_int(count)} simplices...", flush=True)
        if filtration > alpha_sq:
            continue
        if len(simplex) == 2:
            u, v = int(simplex[0]), int(simplex[1])
            if u > v:
                u, v = v, u
            edges.add((u, v))
        elif len(simplex) == 4:
            tet = [int(simplex[0]), int(simplex[1]), int(simplex[2]), int(simplex[3])]
            tetrahedra.append(tet)
            volumes.append(_tetra_volume(pts_xyz[np.asarray(tet, dtype=np.int64)]))
    return edges, tetrahedra, volumes


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--catalog-path", type=Path, required=True, help="Input FITS catalog (e.g. bgs_maglim_*.fits).")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--chunk-rows", type=int, default=2_000_000, help="Rows per streaming chunk.")
    p.add_argument("--alpha-sq", type=float, default=float("inf"), help="Alpha^2 threshold; inf gives full Delaunay.")
    p.add_argument(
        "--output-prefix",
        type=str,
        default=None,
        help="Abacus-style output prefix. Defaults to `desi_delaunay` if alpha_sq is inf else `desi_alpha`.",
    )
    p.add_argument(
        "--simplex-progress-every",
        type=int,
        default=5_000_000,
        help="Print a progress line every N simplices during filtration iteration (0 disables).",
    )
    p.add_argument(
        "--heartbeat-seconds",
        type=int,
        default=60,
        help="Print heartbeat while Gudhi builds simplex tree (0 disables).",
    )
    p.add_argument(
        "--max-points-per-hemi",
        type=int,
        default=0,
        help="Optional cap per hemisphere (0 disables). Randomly subsamples for smoke tests.",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    t0 = time.perf_counter()
    cat_path = args.catalog_path.expanduser().resolve()
    if not cat_path.exists():
        raise FileNotFoundError(cat_path)
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 90, flush=True)
    print("DESI BGS Gudhi graph builder (hemisphere split)", flush=True)
    print("=" * 90, flush=True)
    _print_kv("HOSTNAME", os.environ.get("HOSTNAME", "unknown"))
    _print_kv("SLURM_JOB_ID", os.environ.get("SLURM_JOB_ID", ""))
    _print_kv("SLURM_CPUS_PER_TASK", os.environ.get("SLURM_CPUS_PER_TASK", ""))
    _print_kv("catalog_path", cat_path)
    _print_kv("out_dir", out_dir)
    _print_kv("chunk_rows", _fmt_int(args.chunk_rows))
    _print_kv("alpha_sq", args.alpha_sq)
    _print_kv("max_points_per_hemi", args.max_points_per_hemi)
    _print_kv("seed", args.seed)
    print("-" * 90, flush=True)

    # Columns in the maglim catalog.
    COL_RA = "TARGET_RA"
    COL_DEC = "TARGET_DEC"
    COL_Z = "Z"

    # Pass 1: stream catalog, compute hemisphere and xyz (written to memmap for reproducibility).
    print("Step 1/2: streaming catalog -> hemisphere flag + xyz memmap", flush=True)
    t_pass1 = time.perf_counter()
    with fitsio.FITS(str(cat_path)) as f:
        h = f[1]
        nrows = int(h.get_nrows())
        _print_kv("rows_total", _fmt_int(nrows))

        xyz_path = out_dir / "points_xyz_mpc_h.npy"
        xyz = np.lib.format.open_memmap(xyz_path, mode="w+", dtype=np.float32, shape=(nrows, 3))
        hemi = np.empty((nrows,), dtype=np.int8)  # 1 north, 0 south

        for start in range(0, nrows, int(args.chunk_rows)):
            stop = min(start + int(args.chunk_rows), nrows)
            tab = h.read(rows=range(start, stop), columns=[COL_RA, COL_DEC, COL_Z])
            ra = np.asarray(tab[COL_RA], dtype=np.float64)
            dec = np.asarray(tab[COL_DEC], dtype=np.float64)
            zz = np.asarray(tab[COL_Z], dtype=np.float64)

            sc = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
            b = sc.galactic.b.deg
            hemi[start:stop] = (b > 0).astype(np.int8)

            dist = cosmo.comoving_distance(zz).value * float(cosmo.h)  # Mpc/h
            r = np.deg2rad(ra)
            d = np.deg2rad(dec)
            xx = dist * np.cos(d) * np.cos(r)
            yy = dist * np.cos(d) * np.sin(r)
            zz_c = dist * np.sin(d)
            xyz[start:stop, 0] = xx.astype(np.float32)
            xyz[start:stop, 1] = yy.astype(np.float32)
            xyz[start:stop, 2] = zz_c.astype(np.float32)

            if start == 0 or stop == nrows or (start // int(args.chunk_rows) + 1) % 5 == 0:
                frac = 100.0 * float(stop) / float(nrows)
                print(f"  progress: {frac:6.2f}%  ({_fmt_int(stop)}/{_fmt_int(nrows)})", flush=True)

    hemi_path = out_dir / "hemisphere_flag.npy"
    np.save(hemi_path, hemi)
    _print_kv("wrote_xyz", xyz_path)
    _print_kv("wrote_hemi_flag", hemi_path)
    _print_kv("step1_elapsed", _fmt_sec(time.perf_counter() - t_pass1))
    print("-" * 90, flush=True)

    # Build per-hemisphere edges (Abacus-style).
    print("Step 2/2: Gudhi AlphaComplex per hemisphere -> merged edges", flush=True)
    t_step2 = time.perf_counter()
    try:
        import gudhi  # type: ignore
    except ImportError as exc:
        raise ImportError("Gudhi is required: install `gudhi` in cosmic_env.") from exc

    rng = np.random.default_rng(int(args.seed))
    all_edges: set[tuple[int, int]] = set()
    all_tets: list[list[int]] = []
    all_vols: list[float] = []
    per_hemi = {}
    prefix = args.output_prefix
    if prefix is None:
        prefix = "desi_delaunay" if np.isinf(float(args.alpha_sq)) else "desi_alpha"
    for hemi_flag, hemi_name in [(1, "N"), (0, "S")]:
        t_hemi0 = time.perf_counter()
        idx = np.nonzero(hemi == hemi_flag)[0].astype(np.int64)
        print(f"Hemisphere {hemi_name}: N={_fmt_int(idx.size)}", flush=True)
        if args.max_points_per_hemi and idx.size > int(args.max_points_per_hemi):
            idx = np.sort(rng.choice(idx, size=int(args.max_points_per_hemi), replace=False))
            print(f"  subsample -> {hemi_name}: using {idx.size:,}", flush=True)

        pts = np.asarray(xyz[idx], dtype=np.float64)
        print(f"  building AlphaComplex: points shape={pts.shape}", flush=True)
        t_build0 = time.perf_counter()
        ac = gudhi.AlphaComplex(points=pts)
        st = None
        hb_stop = threading.Event()

        def _heartbeat() -> None:
            if hb_stop.is_set():
                return
            elapsed = time.perf_counter() - t_build0
            print(f"    [heartbeat] still building simplex tree... elapsed={_fmt_sec(elapsed)}", flush=True)
            t = threading.Timer(float(args.heartbeat_seconds), _heartbeat)
            t.daemon = True
            t.start()

        if int(args.heartbeat_seconds) > 0:
            t = threading.Timer(float(args.heartbeat_seconds), _heartbeat)
            t.daemon = True
            t.start()
        st = ac.create_simplex_tree()
        hb_stop.set()
        t_build1 = time.perf_counter()
        # Abacus-parity: always export tetrahedra (3-simplices) + volumes.
        edges_local, tets_local, vols_local = _extract_edges_and_tets(
            st,
            pts,
            float(args.alpha_sq),
            progress_every=int(args.simplex_progress_every),
        )
        t_edges1 = time.perf_counter()
        print(f"  edges_local={_fmt_int(len(edges_local))}", flush=True)
        print(f"  build_elapsed={_fmt_sec(t_build1 - t_build0)}  extract_elapsed={_fmt_sec(t_edges1 - t_build1)}", flush=True)

        # Map local indices back to global ids and accumulate.
        t_map0 = time.perf_counter()
        edges_global = np.empty((len(edges_local), 2), dtype=np.int64)
        for i, (lu, lv) in enumerate(sorted(edges_local)):
            edges_global[i, 0] = idx[lu]
            edges_global[i, 1] = idx[lv]
        edges_global = edges_global.astype(np.int32)
        per_hemi[hemi_name] = int(edges_global.shape[0])
        for eu, ev in edges_global:
            uu = int(eu)
            vv = int(ev)
            if uu > vv:
                uu, vv = vv, uu
            all_edges.add((uu, vv))

        if tets_local:
            for tet in tets_local:
                all_tets.append([int(idx[tet[0]]), int(idx[tet[1]]), int(idx[tet[2]]), int(idx[tet[3]])])
            all_vols.extend(vols_local)
        t_map1 = time.perf_counter()
        print(f"  map+merge_elapsed={_fmt_sec(t_map1 - t_map0)}  hemi_elapsed={_fmt_sec(t_map1 - t_hemi0)}", flush=True)

    edges_arr = np.asarray(sorted(all_edges), dtype=np.int32)
    # Legacy output name (kept for compatibility with existing downstream scripts).
    edges_path_legacy = out_dir / "edges_combined_idx.npy"
    np.save(edges_path_legacy, edges_arr)
    _print_kv("wrote_edges", edges_path_legacy)
    _print_kv("edges_shape", tuple(int(x) for x in edges_arr.shape))
    if edges_arr.size > 0:
        _print_kv("edge_index_range", f"[{int(edges_arr.min())}, {int(edges_arr.max())}]")

    # Legacy names + Abacus-style names.
    tets_arr = np.asarray(all_tets, dtype=np.int32) if all_tets else np.empty((0, 4), dtype=np.int32)
    vols_arr = np.asarray(all_vols, dtype=np.float64) if all_vols else np.empty((0,), dtype=np.float64)
    tet_idx_path_legacy = out_dir / "tetrahedra_idx.npy"
    tet_vol_path_legacy = out_dir / "tetrahedra_volumes.npy"
    np.save(tet_idx_path_legacy, tets_arr)
    np.save(tet_vol_path_legacy, vols_arr)
    _print_kv("wrote_tetrahedra_idx", tet_idx_path_legacy)
    _print_kv("wrote_tetrahedra_volumes", tet_vol_path_legacy)
    _print_kv("tetrahedra_shape", tuple(int(x) for x in tets_arr.shape))

    edges_path_abacus = out_dir / f"{prefix}_edges_combined_idx.npy"
    tet_idx_path_abacus = out_dir / f"{prefix}_tetrahedra_idx.npy"
    tet_vol_path_abacus = out_dir / f"{prefix}_tetrahedra_volumes.npy"
    points_xyz_abacus = out_dir / f"{prefix}_points_xyz.npy"
    meta_path_abacus = out_dir / f"{prefix}_metadata.json"
    np.save(edges_path_abacus, edges_arr)
    np.save(tet_idx_path_abacus, tets_arr)
    np.save(tet_vol_path_abacus, vols_arr)
    # Abacus writes float64 points; here we only have xyz.
    np.save(points_xyz_abacus, np.asarray(xyz, dtype=np.float64))
    _print_kv("wrote_abacus_edges", edges_path_abacus)
    _print_kv("wrote_abacus_tetra_idx", tet_idx_path_abacus)
    _print_kv("wrote_abacus_tetra_vol", tet_vol_path_abacus)
    _print_kv("wrote_abacus_points_xyz", points_xyz_abacus)

    meta = {
        "catalog_path": str(cat_path),
        "nrows": int(nrows),
        "alpha_sq": float(args.alpha_sq),
        "hemisphere_split": True,
        "max_points_per_hemi": int(args.max_points_per_hemi),
        "seed": int(args.seed),
        "points_xyz_path": str(xyz_path),
        "hemisphere_flag_path": str(hemi_path),
        "edges_path": str(edges_path_legacy),
        "edges_per_hemi": per_hemi,
        "total_edges_written": int(edges_arr.shape[0]),
        "tetrahedra_idx_path": str(tet_idx_path_legacy),
        "tetrahedra_volumes_path": str(tet_vol_path_legacy),
        "total_tetrahedra_written": int(tets_arr.shape[0]),
        "abacus_style": {
            "prefix": prefix,
            "edges": str(edges_path_abacus),
            "tetrahedra_idx": str(tet_idx_path_abacus),
            "tetrahedra_volumes": str(tet_vol_path_abacus),
            "points_xyz": str(points_xyz_abacus),
            "metadata": str(meta_path_abacus),
        },
    }
    meta_path = out_dir / "graph_metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _print_kv("wrote_metadata", meta_path)

    # Also write Abacus-style metadata file with the same top-level keys that Abacus emits.
    meta_abacus = {
        "prefix": prefix,
        "mode": "delaunay" if np.isinf(float(args.alpha_sq)) else "alpha",
        "alpha_sq": None if np.isinf(float(args.alpha_sq)) else float(args.alpha_sq),
        "split_hemispheres": True,
        "source": "catalog",
        "source_path": str(cat_path),
        "n_points": int(nrows),
        "n_point_columns": 3,
        "n_edges": int(edges_arr.shape[0]),
        "n_tetrahedra": int(tets_arr.shape[0]),
        "catalog_filters": None,
        "files": {
            "points_xyz": points_xyz_abacus.name,
            "edges": edges_path_abacus.name,
            "tetrahedra_idx": tet_idx_path_abacus.name,
            "tetrahedra_volumes": tet_vol_path_abacus.name,
        },
    }
    meta_path_abacus.write_text(json.dumps(meta_abacus, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _print_kv("wrote_abacus_metadata", meta_path_abacus)
    _print_kv("step2_elapsed", _fmt_sec(time.perf_counter() - t_step2))
    _print_kv("total_elapsed", _fmt_sec(time.perf_counter() - t0))
    print("=" * 90, flush=True)


if __name__ == "__main__":
    main()

