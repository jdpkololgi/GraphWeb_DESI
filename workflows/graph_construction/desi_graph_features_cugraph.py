#!/usr/bin/env python3
"""
Compute DESI graph features (Abacus-style) from Gudhi-built hemisphere graph artifacts.

This mirrors the intent/column names of:
  `TNG/Illustris/workflows/abacus_tweb/abacus_graph_features_cugraph.py`

Node features (7):
  - Degree            : weighted degree = sum(edge_length) per node
  - Clustering        : local clustering coefficient from cuGraph triangle counts
  - Density           : tetrahedral density (matches Abacus definition)
  - Neigh Density     : mean neighbor Density
  - I_eig1/2/3        : eigenvalues of neighbor-position covariance (inertia proxy)

Edge features (5):
  - edge_length, x_dir, y_dir, z_dir, density_contrast

Inputs are discovered from the DESI graph builder metadata JSON produced by:
  `workflows/graph_construction/build_desi_bgs_gudhi_graph.py`

Notes / differences vs Abacus:
  - DESI now stores tetrahedra + volumes (Gudhi 3-simplices), so Density is computed with
    the exact same tetrahedral definition as Abacus:
        node_tetra_count / node_tetra_volume / weighted_degree

  - cuGraph triangle counting + exporting full edge feature tables can be expensive for
    ~1e8 edges. Use output caps if needed.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


NODE_COLS = ["Degree", "Clustering", "Density", "Neigh Density", "I_eig1", "I_eig2", "I_eig3"]
EDGE_COLS = ["edge_length", "x_dir", "y_dir", "z_dir", "density_contrast"]


def _fmt_int(x: int) -> str:
    return f"{int(x):,}"


def _fmt_sec(dt: float) -> str:
    if dt < 120:
        return f"{dt:.1f}s"
    if dt < 7200:
        return f"{dt/60.0:.1f}m"
    return f"{dt/3600.0:.2f}h"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata-path", type=Path, required=True, help="Path to DESI graph_metadata.json.")
    p.add_argument("--out-dir", type=Path, default=None, help="Output directory (default: metadata dir).")
    p.add_argument("--out-prefix", type=str, default="desi_delaunay_cugraph", help="Output prefix.")
    p.add_argument("--batch-size", type=int, default=2_000_000, help="Batch size for edge computations.")
    p.add_argument("--heartbeat-seconds", type=int, default=60, help="Heartbeat print interval (0 disables).")

    p.add_argument(
        "--skip-clustering",
        action="store_true",
        help="Skip cuGraph clustering (set to 0). Useful if cuGraph unavailable or too large.",
    )
    p.add_argument(
        "--skip-inertia",
        action="store_true",
        help="Skip inertia eigs (set to 0). Useful for quick runs.",
    )
    p.add_argument(
        "--write-parquet",
        action="store_true",
        help="Also write node/edge parquet tables (very large for full DESI).",
    )
    p.add_argument(
        "--max-edges-export",
        type=int,
        default=0,
        help="Optional cap on edges exported to parquet/npz edge_attr (0 disables).",
    )
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--debug-inertia-parity",
        action="store_true",
        help=(
            "Run a small parity check comparing inertia eigenvalues computed via "
            "(a) Abacus-style per-node neighbor loop vs (b) vectorized edge-accumulation. "
            "This runs on a random induced subgraph and exits after printing diffs."
        ),
    )
    p.add_argument("--debug-n-nodes", type=int, default=20_000, help="Nodes sampled for inertia debug subgraph.")
    p.add_argument("--debug-n-edges", type=int, default=200_000, help="Edges sampled for inertia debug subgraph.")
    return p.parse_args()


def _load_metadata(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _edge_lengths(points_xyz: np.ndarray, edges: np.ndarray, batch_size: int) -> np.ndarray:
    out = np.empty((edges.shape[0],), dtype=np.float32)
    for i in range(0, edges.shape[0], batch_size):
        b = edges[i : i + batch_size]
        diffs = points_xyz[b[:, 0]] - points_xyz[b[:, 1]]
        out[i : i + b.shape[0]] = np.linalg.norm(diffs, axis=1).astype(np.float32)
    return out


def _degree_unweighted(n_nodes: int, edges: np.ndarray) -> np.ndarray:
    deg = np.zeros((n_nodes,), dtype=np.int32)
    np.add.at(deg, edges[:, 0], 1)
    np.add.at(deg, edges[:, 1], 1)
    return deg


def _degree_weighted(n_nodes: int, edges: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    deg = np.zeros((n_nodes,), dtype=np.float32)
    np.add.at(deg, edges[:, 0], lengths)
    np.add.at(deg, edges[:, 1], lengths)
    return deg


def _tetra_density(n_nodes: int, tetrahedra: np.ndarray, volumes: np.ndarray, degree_w: np.ndarray) -> np.ndarray:
    if tetrahedra.size == 0:
        return np.zeros((n_nodes,), dtype=np.float32)

    node_indices = tetrahedra.reshape(-1)
    repeated_vol = np.repeat(volumes.astype(np.float64, copy=False), 4)

    node_tetra_count = np.bincount(node_indices, minlength=n_nodes).astype(np.float64)
    node_tetra_vol = np.bincount(node_indices, weights=repeated_vol, minlength=n_nodes).astype(np.float64)

    deg_safe = np.where(degree_w > 0, degree_w.astype(np.float64), 1.0)
    dens = np.zeros((n_nodes,), dtype=np.float64)
    m = node_tetra_vol > 0
    dens[m] = (node_tetra_count[m] / node_tetra_vol[m]) / deg_safe[m]
    return dens.astype(np.float32)


def _neighbor_density_mean(n_nodes: int, edges: np.ndarray, density: np.ndarray) -> np.ndarray:
    neigh_sum = np.zeros((n_nodes,), dtype=np.float32)
    neigh_cnt = np.zeros((n_nodes,), dtype=np.int32)
    np.add.at(neigh_sum, edges[:, 0], density[edges[:, 1]])
    np.add.at(neigh_sum, edges[:, 1], density[edges[:, 0]])
    np.add.at(neigh_cnt, edges[:, 0], 1)
    np.add.at(neigh_cnt, edges[:, 1], 1)
    out = np.zeros((n_nodes,), dtype=np.float32)
    m = neigh_cnt > 0
    out[m] = neigh_sum[m] / neigh_cnt[m].astype(np.float32)
    return out


def _cugraph_clustering(n_nodes: int, edges: np.ndarray, degree_u: np.ndarray) -> np.ndarray:
    import cudf  # type: ignore
    import cugraph  # type: ignore

    df = cudf.DataFrame({"src": edges[:, 0].astype(np.int32), "dst": edges[:, 1].astype(np.int32)})
    g = cugraph.Graph(directed=False)
    g.from_cudf_edgelist(df, source="src", destination="dst", renumber=False)

    tc = cugraph.triangle_count(g)
    tri_counts = tc["counts"].to_numpy()
    tri_vertices = tc["vertex"].to_numpy()
    out = np.zeros((n_nodes,), dtype=np.float32)
    out[tri_vertices] = tri_counts.astype(np.float32)

    k = degree_u.astype(np.float32)
    denom = k * (k - 1.0)
    valid = denom > 0
    out[valid] = (2.0 * out[valid]) / denom[valid]
    out[~valid] = 0.0
    return out


def _inertia_eigs_from_edges(points_xyz: np.ndarray, edges: np.ndarray, degree_u: np.ndarray, batch_size: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Vectorized inertia proxy:
      For each node i, take neighbor positions {x_j} and compute covariance about their mean.
    This avoids a Python loop over nodes.
    """
    n = points_xyz.shape[0]
    deg = degree_u.astype(np.int32)
    # We will treat graph as directed both ways.
    src = np.concatenate([edges[:, 0], edges[:, 1]]).astype(np.int64)
    dst = np.concatenate([edges[:, 1], edges[:, 0]]).astype(np.int64)

    # Abacus parity: compute covariance about neighbor mean in float64
    # (Abacus uses float64 points and per-node cov = (rel^T rel)/k).
    sum_pos = np.zeros((n, 3), dtype=np.float64)
    sum_xx = np.zeros((n, 3, 3), dtype=np.float64)

    for i in range(0, src.size, batch_size):
        s = src[i : i + batch_size]
        d = dst[i : i + batch_size]
        p = points_xyz[d].astype(np.float64, copy=False)
        np.add.at(sum_pos, s, p)
        # outer products for each edge endpoint
        # (B,3,3) = (B,3,1) * (B,1,3)
        outer = p[:, :, None] * p[:, None, :]
        np.add.at(sum_xx, s, outer)

    deg_safe = np.where(deg > 0, deg, 1).astype(np.float64)
    mu = sum_pos / deg_safe[:, None]
    exx = sum_xx / deg_safe[:, None, None]
    cov = exx - (mu[:, :, None] * mu[:, None, :])
    # numerical safety
    cov = 0.5 * (cov + np.transpose(cov, (0, 2, 1)))

    # Abacus behavior: if node has <3 neighbors, leave inertia eigs at 0.
    cov[deg < 3] = 0.0

    # eigvalsh in chunks to limit peak memory
    out = np.zeros((n, 3), dtype=np.float32)
    chunk = 1_000_000
    for j in range(0, n, chunk):
        c = cov[j : j + chunk]
        vals = np.linalg.eigvalsh(c).astype(np.float32)
        vals[vals < 0] = 0.0
        out[j : j + vals.shape[0]] = vals
        if j > 0 and j % (5_000_000) == 0:
            print(f"  inertia eigs processed {_fmt_int(j)}/{_fmt_int(n)} nodes", flush=True)

    return out[:, 0], out[:, 1], out[:, 2]


def _inertia_eigs_reference_loop(points_xyz: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Abacus-style reference:
      - build adjacency
      - for each node, compute covariance about neighbor mean using neighbor positions
      - if <3 neighbors => zeros
    """
    from scipy.sparse import csr_matrix

    n_nodes = int(points_xyz.shape[0])
    row = np.concatenate([edges[:, 0], edges[:, 1]]).astype(np.int64)
    col = np.concatenate([edges[:, 1], edges[:, 0]]).astype(np.int64)
    data = np.ones((row.shape[0],), dtype=np.int8)
    adj = csr_matrix((data, (row, col)), shape=(n_nodes, n_nodes))

    eig = np.zeros((n_nodes, 3), dtype=np.float32)
    indptr = adj.indptr
    indices = adj.indices
    pts64 = np.asarray(points_xyz, dtype=np.float64)

    for node in range(n_nodes):
        nbr = indices[indptr[node] : indptr[node + 1]]
        if nbr.size < 3:
            continue
        nbr_pos = pts64[nbr, :3]
        center = nbr_pos.mean(axis=0)
        rel = nbr_pos - center
        cov = (rel.T @ rel) / float(nbr.size)
        vals = np.linalg.eigvalsh(cov).astype(np.float32)
        vals[vals < 0] = 0.0
        eig[node] = vals
    return eig[:, 0], eig[:, 1], eig[:, 2]


def _run_inertia_parity_debug(
    *,
    xyz: np.ndarray,
    edges_full: np.ndarray,
    seed: int,
    n_nodes_sample: int,
    n_edges_sample: int,
    batch_size: int,
) -> None:
    """
    Sample a small induced subgraph and compare inertia eigs computed by:
      - reference loop (Abacus-style)
      - vectorized edge accumulation
    """
    rng = np.random.default_rng(int(seed))
    n_total_nodes = int(xyz.shape[0])
    n_total_edges = int(edges_full.shape[0])

    n_nodes_sample = min(int(n_nodes_sample), n_total_nodes)
    n_edges_sample = min(int(n_edges_sample), n_total_edges)

    node_sel = np.sort(rng.choice(n_total_nodes, size=n_nodes_sample, replace=False))
    node_to_local = {int(g): int(i) for i, g in enumerate(node_sel)}

    # Sample edges until we get enough that lie fully inside the selected node set.
    # This keeps the debug fast and ensures both methods operate on the same graph.
    want = n_edges_sample
    # Use a set to avoid duplicate edges in the sampled subgraph.
    got_set: set[tuple[int, int]] = set()
    trials = 0
    max_trials = 10 * want if want > 0 else 0
    while len(got_set) < want and trials < max_trials:
        batch = rng.integers(
            0,
            n_total_edges,
            size=min(50_000, max(1, want - len(got_set))),
            endpoint=False,
        )
        e = np.asarray(edges_full[batch], dtype=np.int64)
        for u, v in e:
            lu = node_to_local.get(int(u))
            if lu is None:
                continue
            lv = node_to_local.get(int(v))
            if lv is None:
                continue
            if lu > lv:
                lu, lv = lv, lu
            got_set.add((lu, lv))
            if len(got_set) >= want:
                break
        trials += int(batch.size)

    edges = np.asarray(sorted(got_set), dtype=np.int64)
    if edges.size == 0:
        print("Inertia debug: sampled subgraph had 0 internal edges; increase --debug-n-edges or --debug-n-nodes.")
        return

    # Ensure undirected list is consistent (both methods treat as undirected via doubling).
    edges = edges.astype(np.int64, copy=False)
    pts = np.asarray(xyz[node_sel], dtype=np.float64)
    deg_u = _degree_unweighted(int(pts.shape[0]), edges).astype(np.int32)

    print("-" * 90, flush=True)
    print(
        f"Inertia parity debug subgraph: nodes={pts.shape[0]:,} edges={edges.shape[0]:,} "
        f"(sampled from N={n_total_nodes:,} E={n_total_edges:,})",
        flush=True,
    )

    t0 = time.perf_counter()
    r1, r2, r3 = _inertia_eigs_reference_loop(pts, edges)
    t1 = time.perf_counter()
    v1, v2, v3 = _inertia_eigs_from_edges(pts, edges, deg_u, batch_size=batch_size)
    t2 = time.perf_counter()
    print(f"Reference loop elapsed: {_fmt_sec(t1 - t0)}", flush=True)
    print(f"Vectorized elapsed:    {_fmt_sec(t2 - t1)}", flush=True)

    ref = np.stack([r1, r2, r3], axis=1).astype(np.float64)
    vec = np.stack([v1, v2, v3], axis=1).astype(np.float64)
    diff = np.abs(ref - vec)
    max_diff = float(np.max(diff))
    p99 = float(np.percentile(diff, 99))
    p50 = float(np.percentile(diff, 50))
    print(f"abs diff stats: max={max_diff:.3e}  p99={p99:.3e}  p50={p50:.3e}", flush=True)
    # Also check deg<3 masking matches
    mask = deg_u < 3
    if np.any(mask):
        nonzero_ref = float(np.max(np.abs(ref[mask])))
        nonzero_vec = float(np.max(np.abs(vec[mask])))
        print(f"deg<3 check: max|ref|={nonzero_ref:.3e}  max|vec|={nonzero_vec:.3e}", flush=True)
    print("-" * 90, flush=True)


def main() -> None:
    args = parse_args()
    t0 = time.perf_counter()

    meta_path = args.metadata_path.expanduser().resolve()
    meta = _load_metadata(meta_path)
    out_dir = (args.out_dir.expanduser().resolve() if args.out_dir else meta_path.parent.resolve())
    out_dir.mkdir(parents=True, exist_ok=True)

    points_path = Path(meta["points_xyz_path"]).expanduser().resolve()
    edges_path = Path(meta["edges_path"]).expanduser().resolve()
    tetra_path_raw = meta.get("tetrahedra_idx_path")
    vol_path_raw = meta.get("tetrahedra_volumes_path")
    if tetra_path_raw is None or vol_path_raw is None:
        raise ValueError(
            "Metadata is missing tetrahedra paths. Rebuild graph with updated builder so "
            "tetrahedra_idx/tetrahedra_volumes are exported."
        )
    tetra_path = Path(tetra_path_raw).expanduser().resolve()
    vol_path = Path(vol_path_raw).expanduser().resolve()
    if not points_path.exists():
        raise FileNotFoundError(points_path)
    if not edges_path.exists():
        raise FileNotFoundError(edges_path)
    if not tetra_path.exists():
        raise FileNotFoundError(tetra_path)
    if not vol_path.exists():
        raise FileNotFoundError(vol_path)

    print("=" * 90, flush=True)
    print("DESI graph features (cuGraph-style)", flush=True)
    print("=" * 90, flush=True)
    print("HOSTNAME:", os.environ.get("HOSTNAME", "unknown"), flush=True)
    print("SLURM_JOB_ID:", os.environ.get("SLURM_JOB_ID", ""), flush=True)
    print("metadata_path:", meta_path, flush=True)
    print("points_xyz_path:", points_path, flush=True)
    print("edges_path:", edges_path, flush=True)
    print("tetrahedra_idx_path:", tetra_path, flush=True)
    print("tetrahedra_volumes_path:", vol_path, flush=True)
    print("out_dir:", out_dir, flush=True)
    print("out_prefix:", args.out_prefix, flush=True)
    print("-" * 90, flush=True)

    pts = np.load(points_path, mmap_mode="r")
    edges = np.load(edges_path, mmap_mode="r")
    tetra = np.load(tetra_path, mmap_mode="r")
    vols = np.load(vol_path, mmap_mode="r")
    n_nodes = int(pts.shape[0])
    n_edges = int(edges.shape[0])
    print(f"Loaded points={pts.shape} dtype={pts.dtype}", flush=True)
    print(f"Loaded edges={edges.shape} dtype={edges.dtype}", flush=True)
    print(f"Loaded tetrahedra={tetra.shape} dtype={tetra.dtype}", flush=True)
    print(f"Loaded tetra volumes={vols.shape} dtype={vols.dtype}", flush=True)

    # For Abacus parity: compute ALL node-level metrics on the FULL edge list.
    # Edge capping (if enabled) applies only to exported edge attributes/parquet,
    # never to node metrics.
    edges_full = edges
    export_edges = edges_full
    if args.max_edges_export and int(edges_full.shape[0]) > int(args.max_edges_export):
        rng = np.random.default_rng(int(args.seed))
        sel = rng.choice(int(edges_full.shape[0]), size=int(args.max_edges_export), replace=False)
        export_edges = np.asarray(edges_full[np.sort(sel)], dtype=np.int64)
        print(f"Applied --max-edges-export: using edges={export_edges.shape[0]:,} for edge-attr export", flush=True)
    else:
        # Keep export edges as a materialized int64 array for consistent downstream dtypes.
        export_edges = np.asarray(edges_full, dtype=np.int64)

    # Materialize xyz as float32 for computations.
    xyz = np.asarray(pts[:, :3], dtype=np.float32)

    if args.debug_inertia_parity:
        _run_inertia_parity_debug(
            xyz=xyz,
            edges_full=np.asarray(edges, dtype=np.int64),
            seed=int(args.seed),
            n_nodes_sample=int(args.debug_n_nodes),
            n_edges_sample=int(args.debug_n_edges),
            batch_size=int(args.batch_size),
        )
        return

    print("-" * 90, flush=True)
    print("Computing edge lengths...", flush=True)
    t_len0 = time.perf_counter()
    lengths_full = _edge_lengths(xyz, np.asarray(edges_full, dtype=np.int64), int(args.batch_size))
    print(f"  lengths computed in {_fmt_sec(time.perf_counter() - t_len0)}", flush=True)

    print("Computing degrees...", flush=True)
    t_deg0 = time.perf_counter()
    degree_u = _degree_unweighted(n_nodes, np.asarray(edges_full, dtype=np.int64))
    degree_w = _degree_weighted(n_nodes, np.asarray(edges_full, dtype=np.int64), lengths_full)
    print(f"  degrees computed in {_fmt_sec(time.perf_counter() - t_deg0)}", flush=True)

    print("Computing tetrahedral density (Abacus definition)...", flush=True)
    t_den0 = time.perf_counter()
    density = _tetra_density(n_nodes, np.asarray(tetra, dtype=np.int64), np.asarray(vols, dtype=np.float64), degree_w)
    print(f"  density computed in {_fmt_sec(time.perf_counter() - t_den0)}", flush=True)

    print("Computing neighbor density mean...", flush=True)
    t_nd0 = time.perf_counter()
    neigh_density = _neighbor_density_mean(n_nodes, np.asarray(edges_full, dtype=np.int64), density)
    print(f"  neigh density computed in {_fmt_sec(time.perf_counter() - t_nd0)}", flush=True)

    if args.skip_clustering:
        clustering = np.zeros((n_nodes,), dtype=np.float32)
        print("Skipping clustering (set to 0).", flush=True)
    else:
        print("Computing clustering with cuGraph triangle counts...", flush=True)
        t_cl0 = time.perf_counter()
        clustering = _cugraph_clustering(n_nodes, np.asarray(edges_full, dtype=np.int64), degree_u)
        print(f"  clustering computed in {_fmt_sec(time.perf_counter() - t_cl0)}", flush=True)

    if args.skip_inertia:
        i1 = np.zeros((n_nodes,), dtype=np.float32)
        i2 = np.zeros((n_nodes,), dtype=np.float32)
        i3 = np.zeros((n_nodes,), dtype=np.float32)
        print("Skipping inertia eigs (set to 0).", flush=True)
    else:
        print("Computing inertia eigs (neighbor-position covariance)...", flush=True)
        t_i0 = time.perf_counter()
        i1, i2, i3 = _inertia_eigs_from_edges(xyz, np.asarray(edges_full, dtype=np.int64), degree_u, int(args.batch_size))
        print(f"  inertia eigs computed in {_fmt_sec(time.perf_counter() - t_i0)}", flush=True)

    # Edge attributes for export set (may be full or capped).
    lengths = _edge_lengths(xyz, export_edges, int(args.batch_size))
    vec = xyz[export_edges[:, 1]] - xyz[export_edges[:, 0]]
    denom = np.where(lengths > 0, lengths, 1.0).astype(np.float32)
    unit = (vec / denom[:, None]).astype(np.float32)
    src_d = density[export_edges[:, 0]]
    dst_d = density[export_edges[:, 1]]
    density_contrast = np.zeros_like(lengths, dtype=np.float32)
    mm = src_d > 0
    density_contrast[mm] = dst_d[mm] / src_d[mm]

    # Save compact GNN-ready arrays (npz) + optional torch bundle.
    out_npz = out_dir / f"{args.out_prefix}_gnn_arrays.npz"
    x = np.stack([degree_w, clustering, density, neigh_density, i1, i2, i3], axis=1).astype(np.float32)
    edge_index = export_edges.astype(np.int64).T
    edge_attr = np.stack([lengths, unit[:, 0], unit[:, 1], unit[:, 2], density_contrast], axis=1).astype(np.float32)
    np.savez_compressed(out_npz, x=x, edge_index=edge_index, edge_attr=edge_attr)
    print(f"Saved GNN arrays: {out_npz}", flush=True)

    try:
        import torch  # type: ignore

        pt_path = out_npz.with_suffix(".pt")
        torch.save(
            {
                "x": torch.from_numpy(x),
                "edge_index": torch.from_numpy(edge_index),
                "edge_attr": torch.from_numpy(edge_attr),
            },
            pt_path,
        )
        print(f"Saved PyTorch tensor bundle: {pt_path}", flush=True)
    except Exception:
        pass

    if args.write_parquet:
        import pandas as pd

        node_path = out_dir / f"{args.out_prefix}_node_features.parquet"
        edge_path = out_dir / f"{args.out_prefix}_edge_features.parquet"
        node_df = pd.DataFrame(
            {
                "Node ID": np.arange(n_nodes, dtype=np.int64),
                "Degree": degree_w,
                "Clustering": clustering,
                "Density": density,
                "Neigh Density": neigh_density,
                "I_eig1": i1,
                "I_eig2": i2,
                "I_eig3": i3,
            }
        )
        edge_df = pd.DataFrame(
            {
                "src": export_edges[:, 0].astype(np.int64),
                "dst": export_edges[:, 1].astype(np.int64),
                "edge_length": lengths,
                "x_dir": unit[:, 0],
                "y_dir": unit[:, 1],
                "z_dir": unit[:, 2],
                # Abacus exports both names; keep both for parity.
                "density_contrast_along_edge": density_contrast,
                "density_contrast": density_contrast,
            }
        )
        node_df.to_parquet(node_path, index=False)
        edge_df.to_parquet(edge_path, index=False)
        print(f"Saved node features parquet: {node_path}", flush=True)
        print(f"Saved edge features parquet: {edge_path}", flush=True)

    meta_out = out_dir / f"{args.out_prefix}_gnn_metadata.json"
    payload = {
        "input_metadata_path": str(meta_path),
        "points_xyz_path": str(points_path),
        "edges_path": str(edges_path),
        "tetrahedra_idx_path": str(tetra_path),
        "tetrahedra_volumes_path": str(vol_path),
        "n_points": int(n_nodes),
        "n_edges_total": int(n_edges),
        "n_edges_exported": int(export_edges.shape[0]),
        "node_feature_columns": NODE_COLS,
        "edge_feature_columns": EDGE_COLS,
        "density_definition": "density = (node_tetra_count / node_tetra_volume) / weighted_degree",
        "outputs": {
            "gnn_arrays_npz": str(out_npz),
        },
    }
    with meta_out.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    print(f"Saved metadata: {meta_out}", flush=True)
    print(f"Total elapsed: {_fmt_sec(time.perf_counter() - t0)}", flush=True)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()

