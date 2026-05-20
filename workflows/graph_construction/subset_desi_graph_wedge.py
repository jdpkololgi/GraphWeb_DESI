#!/usr/bin/env python3
"""Subset an existing DESI full graph to a wedge (RA/DEC/Z), Abacus-style.

This script is intentionally shaped like:
`TNG/Illustris/workflows/abacus_tweb/subset_abacus_graph_wedge_for_sbi.py`.

It **cuts** (induces) a subgraph from an already-built full DESI graph and writes
the same core artifact names Abacus uses:

  <out-prefix>_points.npy
  <out-prefix>_points_xyz.npy
  <out-prefix>_edges_combined_idx.npy
  <out-prefix>_tetrahedra_idx.npy
  <out-prefix>_tetrahedra_volumes.npy
  <out-prefix>_metadata.json
  <out-prefix>_global_node_ids.npy   (mapping local->original global node id)

Additionally, since DESI already has full-graph cugraph metrics, this script can
subset those metrics (node+edge) using the wedge indices instead of recomputing.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import fitsio


def _fmt_int(x: int) -> str:
    return f"{int(x):,}"


def _fmt_sec(dt: float) -> str:
    if dt < 120:
        return f"{dt:.1f}s"
    if dt < 7200:
        return f"{dt/60.0:.1f}m"
    return f"{dt/3600.0:.2f}h"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _in_range(arr: np.ndarray, lo: float | None, hi: float | None) -> np.ndarray:
    m = np.ones(arr.shape[0], dtype=bool)
    if lo is not None:
        m &= arr >= float(lo)
    if hi is not None:
        m &= arr <= float(hi)
    return m


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--graph-metadata",
        type=Path,
        required=True,
        help="Path to <prefix>_metadata.json produced by build_desi_bgs_gudhi_graph.py (Abacus-style).",
    )
    p.add_argument(
        "--catalog-path",
        type=Path,
        required=True,
        help="Full DESI FITS catalog with TARGET_RA/TARGET_DEC/Z (row order must match full graph).",
    )
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--out-prefix", type=str, required=True)
    p.add_argument("--chunk-rows", type=int, default=2_000_000, help="Catalog rows per chunk.")
    p.add_argument("--edge-chunk", type=int, default=5_000_000, help="Edges per chunk.")
    p.add_argument("--tet-chunk", type=int, default=2_000_000, help="Tetrahedra per chunk.")

    p.add_argument("--ra-min", type=float, required=True)
    p.add_argument("--ra-max", type=float, required=True)
    p.add_argument("--dec-min", type=float, required=True)
    p.add_argument("--dec-max", type=float, required=True)
    p.add_argument("--z-min", type=float, required=True)
    p.add_argument("--z-max", type=float, required=True)

    p.add_argument(
        "--max-nodes",
        type=int,
        default=0,
        help="Optional cap on selected nodes for smoke tests (0 disables). Keeps first N by catalog order.",
    )
    p.add_argument(
        "--parent-gnn-arrays",
        type=Path,
        default=None,
        help="Optional full-graph <...>_gnn_arrays.npz to subset into the wedge.",
    )
    p.add_argument(
        "--parent-gnn-metadata",
        type=Path,
        default=None,
        help="Optional full-graph <...>_gnn_metadata.json to carry through + update paths.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    t0 = time.perf_counter()

    cat_path = args.catalog_path.expanduser().resolve()
    meta_path = args.graph_metadata.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    if not cat_path.exists():
        raise FileNotFoundError(cat_path)
    if not meta_path.exists():
        raise FileNotFoundError(meta_path)

    meta = _load_json(meta_path)
    base_dir = meta_path.parent
    files = meta.get("files", {})
    # DESI Abacus-style metadata always provides points_xyz/edges/tets/vols in files dict.
    points_xyz_path = base_dir / files.get("points_xyz", f"{meta.get('prefix')}_points_xyz.npy")
    edges_path = base_dir / files.get("edges", f"{meta.get('prefix')}_edges_combined_idx.npy")
    tetra_path = base_dir / files.get("tetrahedra_idx", f"{meta.get('prefix')}_tetrahedra_idx.npy")
    vol_path = base_dir / files.get("tetrahedra_volumes", f"{meta.get('prefix')}_tetrahedra_volumes.npy")

    for pth in (points_xyz_path, edges_path, tetra_path, vol_path):
        if not pth.exists():
            raise FileNotFoundError(pth)

    print("=" * 90, flush=True)
    print("DESI graph wedge subsetting (cut existing full graph)", flush=True)
    print("=" * 90, flush=True)
    print("catalog_path:", cat_path, flush=True)
    print("graph_metadata_path:", meta_path, flush=True)
    print("out_dir:", out_dir, flush=True)
    print("out_prefix:", args.out_prefix, flush=True)
    print(
        "cuts:",
        f"RA[{args.ra_min},{args.ra_max}] Dec[{args.dec_min},{args.dec_max}] z[{args.z_min},{args.z_max}]",
        flush=True,
    )
    print("-" * 90, flush=True)

    # Step 1: select catalog row indices for the wedge.
    COL_RA = "TARGET_RA"
    COL_DEC = "TARGET_DEC"
    COL_Z = "Z"

    print("Step 1/3: selecting node ids (catalog row indices) by RA/Dec/z", flush=True)
    t_sel0 = time.perf_counter()
    sel_rows: list[np.ndarray] = []
    sel_ra: list[np.ndarray] = []
    sel_dec: list[np.ndarray] = []
    sel_z: list[np.ndarray] = []
    with fitsio.FITS(str(cat_path)) as f:
        h = f[1]
        nrows = int(h.get_nrows())
        for start in range(0, nrows, int(args.chunk_rows)):
            stop = min(start + int(args.chunk_rows), nrows)
            tab = h.read(rows=range(start, stop), columns=[COL_RA, COL_DEC, COL_Z])
            ra = np.asarray(tab[COL_RA], dtype=np.float64)
            dec = np.asarray(tab[COL_DEC], dtype=np.float64)
            zz = np.asarray(tab[COL_Z], dtype=np.float64)

            m = _in_range(ra, args.ra_min, args.ra_max)
            m &= _in_range(dec, args.dec_min, args.dec_max)
            m &= _in_range(zz, args.z_min, args.z_max)
            if np.any(m):
                idx = (np.nonzero(m)[0].astype(np.int64) + int(start))
                sel_rows.append(idx)
                sel_ra.append(ra[m].astype(np.float32))
                sel_dec.append(dec[m].astype(np.float32))
                sel_z.append(zz[m].astype(np.float32))

    if not sel_rows:
        raise RuntimeError("No galaxies matched wedge cuts.")

    wedge_global_ids = np.concatenate(sel_rows)
    ra_sel = np.concatenate(sel_ra)
    dec_sel = np.concatenate(sel_dec)
    z_sel = np.concatenate(sel_z)

    if args.max_nodes and wedge_global_ids.size > int(args.max_nodes):
        wedge_global_ids = wedge_global_ids[: int(args.max_nodes)]
        ra_sel = ra_sel[: int(args.max_nodes)]
        dec_sel = dec_sel[: int(args.max_nodes)]
        z_sel = z_sel[: int(args.max_nodes)]

    global_ids_out = out_dir / f"{args.out_prefix}_global_node_ids.npy"
    np.save(global_ids_out, wedge_global_ids.astype(np.int64, copy=False))
    # Optional minimal catalog fields for quick plotting/debug (not FITS).
    np.savez_compressed(out_dir / f"{args.out_prefix}_wedge_catalog_minimal.npz", ra=ra_sel, dec=dec_sel, z=z_sel)
    print(f"  selected nodes: {_fmt_int(wedge_global_ids.size)}", flush=True)
    print(f"  step1_elapsed: {_fmt_sec(time.perf_counter() - t_sel0)}", flush=True)
    print("-" * 90, flush=True)

    # Step 2: subset node arrays (xyz) and write Abacus-style points/xyz.
    print("Step 2/3: writing wedge points + xyz arrays", flush=True)
    t_nodes0 = time.perf_counter()
    xyz_full = np.load(points_xyz_path, mmap_mode="r")
    wedge_xyz = np.asarray(xyz_full[wedge_global_ids], dtype=np.float64)
    wedge_points = wedge_xyz  # Abacus has a richer points.npy; for DESI we store xyz only.
    out_points = out_dir / f"{args.out_prefix}_points.npy"
    out_xyz = out_dir / f"{args.out_prefix}_points_xyz.npy"
    np.save(out_points, wedge_points)
    np.save(out_xyz, wedge_xyz)
    print(f"  wrote points: {out_points.name} {wedge_points.shape}", flush=True)
    print(f"  wrote xyz:    {out_xyz.name} {wedge_xyz.shape}", flush=True)
    print(f"  step2_elapsed: {_fmt_sec(time.perf_counter() - t_nodes0)}", flush=True)
    print("-" * 90, flush=True)

    # Step 3: stream-subset edges + tetrahedra and remap indices.
    print("Step 3/3: subsetting edges + tetrahedra and remapping indices", flush=True)
    t_graph0 = time.perf_counter()
    n_full = int(xyz_full.shape[0])
    in_set = np.zeros((n_full,), dtype=bool)
    in_set[wedge_global_ids.astype(np.int64, copy=False)] = True
    new_id = np.full((n_full,), -1, dtype=np.int64)
    new_id[wedge_global_ids.astype(np.int64, copy=False)] = np.arange(wedge_global_ids.size, dtype=np.int64)

    # --- edges (two-pass like Abacus, streaming) ---
    edges_full = np.load(edges_path, mmap_mode="r")
    if edges_full.ndim != 2 or edges_full.shape[1] != 2:
        raise ValueError(f"edges must be (E,2); got {edges_full.shape}")
    n_keep = 0
    for i in range(0, edges_full.shape[0], int(args.edge_chunk)):
        e = edges_full[i : i + int(args.edge_chunk)]
        m = in_set[e[:, 0]] & in_set[e[:, 1]]
        n_keep += int(np.count_nonzero(m))
    out_edges = out_dir / f"{args.out_prefix}_edges_combined_idx.npy"
    edges_out = np.empty((n_keep, 2), dtype=np.int32)
    w = 0
    for i in range(0, edges_full.shape[0], int(args.edge_chunk)):
        e = edges_full[i : i + int(args.edge_chunk)]
        m = in_set[e[:, 0]] & in_set[e[:, 1]]
        if not np.any(m):
            continue
        sel = e[m].astype(np.int64, copy=False)
        mapped = new_id[sel].astype(np.int32, copy=False)
        edges_out[w : w + mapped.shape[0]] = mapped
        w += mapped.shape[0]
    if w != n_keep:
        raise RuntimeError(f"edge write mismatch: wrote {w} != expected {n_keep}")
    np.save(out_edges, edges_out)

    # --- tets (two-pass like Abacus, streaming) ---
    tets_full = np.load(tetra_path, mmap_mode="r")
    vols_full = np.load(vol_path, mmap_mode="r")
    if tets_full.ndim != 2 or tets_full.shape[1] != 4:
        raise ValueError(f"tetrahedra must be (T,4); got {tets_full.shape}")
    if vols_full.ndim != 1 or vols_full.shape[0] != tets_full.shape[0]:
        raise ValueError("tetrahedra volumes must be (T,) and match tetrahedra rows")
    n_keep_t = 0
    for i in range(0, tets_full.shape[0], int(args.tet_chunk)):
        t = tets_full[i : i + int(args.tet_chunk)]
        m = in_set[t[:, 0]] & in_set[t[:, 1]] & in_set[t[:, 2]] & in_set[t[:, 3]]
        n_keep_t += int(np.count_nonzero(m))
    out_tets = out_dir / f"{args.out_prefix}_tetrahedra_idx.npy"
    out_vols = out_dir / f"{args.out_prefix}_tetrahedra_volumes.npy"
    tets_out = np.empty((n_keep_t, 4), dtype=np.int32)
    vols_out = np.empty((n_keep_t,), dtype=np.float64)
    w = 0
    for i in range(0, tets_full.shape[0], int(args.tet_chunk)):
        t = tets_full[i : i + int(args.tet_chunk)]
        m = in_set[t[:, 0]] & in_set[t[:, 1]] & in_set[t[:, 2]] & in_set[t[:, 3]]
        if not np.any(m):
            continue
        sel_t = t[m].astype(np.int64, copy=False)
        mapped_t = new_id[sel_t].astype(np.int32, copy=False)
        vv = np.asarray(vols_full[i : i + int(args.tet_chunk)][m], dtype=np.float64)
        tets_out[w : w + mapped_t.shape[0]] = mapped_t
        vols_out[w : w + vv.shape[0]] = vv
        w += mapped_t.shape[0]
    if w != n_keep_t:
        raise RuntimeError(f"tet write mismatch: wrote {w} != expected {n_keep_t}")
    np.save(out_tets, tets_out)
    np.save(out_vols, vols_out)

    print(f"  wedge nodes: {_fmt_int(wedge_global_ids.size)} / {_fmt_int(n_full)}", flush=True)
    print(f"  wedge edges: {_fmt_int(edges_out.shape[0])}", flush=True)
    print(f"  wedge tetrahedra: {_fmt_int(tets_out.shape[0])}", flush=True)
    print(f"  step3_elapsed: {_fmt_sec(time.perf_counter() - t_graph0)}", flush=True)

    # Write Abacus-style wedge metadata.
    meta_out = out_dir / f"{args.out_prefix}_metadata.json"
    out_meta = {
        "prefix": args.out_prefix,
        "parent_graph_metadata": str(meta_path),
        "mode": meta.get("mode"),
        "alpha_sq": meta.get("alpha_sq"),
        "split_hemispheres": bool(meta.get("split_hemispheres", True)),
        "source": "catalog",
        "source_path": str(cat_path),
        "n_points": int(wedge_points.shape[0]),
        "n_point_columns": int(wedge_points.shape[1]),
        "n_edges": int(edges_out.shape[0]),
        "n_tetrahedra": int(tets_out.shape[0]),
        "catalog_filters": {
            "ra_min": float(args.ra_min),
            "ra_max": float(args.ra_max),
            "dec_min": float(args.dec_min),
            "dec_max": float(args.dec_max),
            "z_min": float(args.z_min),
            "z_max": float(args.z_max),
        },
        "files": {
            "points": out_points.name,
            "points_xyz": out_xyz.name,
            "edges": out_edges.name,
            "tetrahedra_idx": out_tets.name,
            "tetrahedra_volumes": out_vols.name,
            "global_node_ids": global_ids_out.name,
        },
    }
    meta_out.write_text(json.dumps(out_meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"  wrote wedge graph metadata: {meta_out}", flush=True)

    # Optional: subset existing GNN arrays/metadata using wedge indices.
    if args.parent_gnn_arrays is not None:
        npz_in = args.parent_gnn_arrays.expanduser().resolve()
        if not npz_in.exists():
            raise FileNotFoundError(npz_in)
        data = np.load(npz_in)
        x_full = np.asarray(data["x"])
        edge_index_full = np.asarray(data["edge_index"])  # (2,E)
        edge_attr_full = np.asarray(data["edge_attr"])

        x_sub = x_full[wedge_global_ids.astype(np.int64, copy=False)]
        # Keep only edges fully inside wedge, remap to local ids.
        src = edge_index_full[0].astype(np.int64, copy=False)
        dst = edge_index_full[1].astype(np.int64, copy=False)
        m = in_set[src] & in_set[dst]
        src_l = new_id[src[m]].astype(np.int64, copy=False)
        dst_l = new_id[dst[m]].astype(np.int64, copy=False)
        edge_index_sub = np.stack([src_l, dst_l], axis=0)
        edge_attr_sub = edge_attr_full[m]

        out_npz = out_dir / f"{args.out_prefix}_gnn_arrays.npz"
        np.savez_compressed(out_npz, x=x_sub.astype(np.float32), edge_index=edge_index_sub, edge_attr=edge_attr_sub.astype(np.float32))
        print(f"  wrote subset GNN arrays: {out_npz}", flush=True)

        if args.parent_gnn_metadata is not None:
            gmeta_in = args.parent_gnn_metadata.expanduser().resolve()
            if not gmeta_in.exists():
                raise FileNotFoundError(gmeta_in)
            gmeta = _load_json(gmeta_in)
            gmeta_out = out_dir / f"{args.out_prefix}_gnn_metadata.json"
            gmeta["input_parent_gnn_metadata_path"] = str(gmeta_in)
            gmeta["input_wedge_graph_metadata_path"] = str(meta_out)
            gmeta["n_points"] = int(x_sub.shape[0])
            gmeta["n_edges_total"] = int(edge_index_sub.shape[1])
            gmeta["n_edges_exported"] = int(edge_index_sub.shape[1])
            gmeta.setdefault("outputs", {})
            gmeta["outputs"]["gnn_arrays_npz"] = str(out_npz)
            gmeta_out.write_text(json.dumps(gmeta, indent=2, sort_keys=True) + "\n", encoding="utf-8")
            print(f"  wrote subset GNN metadata: {gmeta_out}", flush=True)

    print("-" * 90, flush=True)
    print(f"Total elapsed: {_fmt_sec(time.perf_counter() - t0)}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()

