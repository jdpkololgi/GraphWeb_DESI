#!/usr/bin/env python3
"""Build a Jraph-compatible inference cache for a DESI BGS wedge.

This script is designed to run *inference* with the Abacus-trained Jraph regression models.
To avoid constructing an all-sky graph, it reuses the cached full DESI Delaunay graph and
feature tables produced by the GraphWeb_DESI workflow, then extracts an induced subgraph
for a user-defined wedge in (RA, DEC, Z).

Inputs (defaults are from the GraphWeb_DESI repo):
- networkx graph: GraphWeb_DESI/cache/DESI_delaunay_graph.pt
- features dataframe: GraphWeb_DESI/cache/DESI_delaunay_features.pt
- zcat dataframe: GraphWeb_DESI/cache/DESI_NETWORK_delaunay_zcat.pt

Output pickle contains:
- graph: jraph.GraphsTuple (nodes + edges + senders/receivers)
- node_feature_names: list[str] (7 columns expected by Abacus-trained models)
- galaxy_table: pandas.DataFrame (subset, aligned to graph node order)
- xyz_mpc: np.ndarray [N,3] float64 comoving Mpc (derived from RA/DEC/Z; Abacus parity)
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Allow canonical workflow scripts to resolve repo-root modules after reorganization.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import jraph


NODE_FEATURE_NAMES = (
    "Degree",
    "Clustering",
    "Density",
    "Neigh Density",
    "I_eig1",
    "I_eig2",
    "I_eig3",
)


def _default_repo_root() -> Path:
    return Path(os.environ.get("GRAPHWEB_REPO_ROOT", str(REPO_ROOT)))


def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ra-min", type=float, required=True)
    p.add_argument("--ra-max", type=float, required=True)
    p.add_argument("--dec-min", type=float, required=True)
    p.add_argument("--dec-max", type=float, required=True)
    p.add_argument("--z-min", type=float, required=True)
    p.add_argument("--z-max", type=float, required=True)
    p.add_argument(
        "--repo-root",
        type=Path,
        default=_default_repo_root(),
        help="GraphWeb_DESI repo root (used to locate cache/ by default).",
    )
    p.add_argument(
        "--out-cache",
        type=Path,
        required=True,
        help="Output pickle path for the wedge inference cache.",
    )
    p.add_argument(
        "--use-z-col",
        type=str,
        default="Z",
        help="Column in zcat to use for redshift cuts (default: Z).",
    )
    p.add_argument(
        "--no-bidirectional-edges",
        action="store_true",
        help="Do not duplicate reverse-direction edges (default: duplicate, matching Abacus caches).",
    )
    p.add_argument(
        "--max-nodes",
        type=int,
        default=0,
        help="Optional cap on number of wedge nodes (0 disables). If set, randomly subsamples wedge nodes for smoke tests.",
    )
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


def _radec_z_to_xyz_mpc(ra_deg: np.ndarray, dec_deg: np.ndarray, z: np.ndarray) -> np.ndarray:
    from shared.abacus_gnn_parity import sky_to_xyz_mpc

    x, y, zc = sky_to_xyz_mpc(ra_deg, dec_deg, z)
    return np.stack([x, y, zc], axis=-1)


def main() -> None:
    args = _parse_args()
    repo = args.repo_root.expanduser().resolve()

    graph_path = repo / "cache" / "DESI_delaunay_graph.pt"
    features_path = repo / "cache" / "DESI_delaunay_features.pt"
    zcat_path = repo / "cache" / "DESI_NETWORK_delaunay_zcat.pt"

    for p in (graph_path, features_path, zcat_path):
        if not p.exists():
            raise FileNotFoundError(f"Missing required GraphWeb_DESI cache file: {p}")

    import torch
    import networkx as nx

    print(f"Loading full DESI Delaunay graph: {graph_path}")
    G_full = torch.load(graph_path, weights_only=False)
    if not isinstance(G_full, nx.Graph):
        raise TypeError(f"Expected networkx Graph at {graph_path}, got {type(G_full)}")

    print(f"Loading DESI features: {features_path}")
    feats_full = pd.read_pickle(features_path)
    missing = [c for c in NODE_FEATURE_NAMES if c not in feats_full.columns]
    if missing:
        raise KeyError(f"Missing required node feature columns in {features_path}: {missing}")

    print(f"Loading DESI zcat: {zcat_path}")
    zcat_full = pd.read_pickle(zcat_path)
    for c in ("RA", "DEC", args.use_z_col):
        if c not in zcat_full.columns:
            raise KeyError(f"zcat missing required column {c!r}.")

    ra = np.asarray(zcat_full["RA"], dtype=np.float64)
    dec = np.asarray(zcat_full["DEC"], dtype=np.float64)
    zz = np.asarray(zcat_full[args.use_z_col], dtype=np.float64)
    m = (
        (ra >= args.ra_min)
        & (ra <= args.ra_max)
        & (dec >= args.dec_min)
        & (dec <= args.dec_max)
        & (zz >= args.z_min)
        & (zz <= args.z_max)
        & np.isfinite(ra)
        & np.isfinite(dec)
        & np.isfinite(zz)
    )
    idx = np.nonzero(m)[0].astype(np.int64)
    print(f"Wedge nodes (pre-cap): {idx.size:,} / {len(zcat_full):,}")
    if idx.size == 0:
        raise ValueError("Wedge selection produced zero nodes; check bounds.")

    if args.max_nodes and idx.size > args.max_nodes:
        rng = np.random.default_rng(int(args.seed))
        idx = np.sort(rng.choice(idx, size=int(args.max_nodes), replace=False))
        print(f"Applied --max-nodes cap: using {idx.size:,} nodes.")

    nodes_old = idx.tolist()
    print("Building induced subgraph...")
    G_sub = G_full.subgraph(nodes_old).copy()
    print(f"Subgraph: nodes={G_sub.number_of_nodes():,}, edges={G_sub.number_of_edges():,}")

    old_to_new = {old: i for i, old in enumerate(nodes_old)}
    N = len(nodes_old)

    feats_sub = feats_full.iloc[nodes_old][list(NODE_FEATURE_NAMES)]
    nodes = np.asarray(feats_sub.values, dtype=np.float32)

    gal_sub = zcat_full.iloc[nodes_old].copy()
    gal_sub.insert(0, "NODE_ID_OLD", np.asarray(nodes_old, dtype=np.int64))
    gal_sub.insert(1, "NODE_ID", np.arange(N, dtype=np.int64))

    xyz = _radec_z_to_xyz_mpc(
        gal_sub["RA"].to_numpy(), gal_sub["DEC"].to_numpy(), gal_sub[args.use_z_col].to_numpy()
    )

    senders = []
    receivers = []
    e_feat = []
    dens = np.asarray(feats_sub["Density"].to_numpy(), dtype=np.float64)
    for u_old, v_old, dat in G_sub.edges(data=True):
        u = old_to_new[int(u_old)]
        v = old_to_new[int(v_old)]
        duv = xyz[v] - xyz[u]
        length = float(np.linalg.norm(duv))
        if not np.isfinite(length) or length <= 0:
            continue
        ddir = duv / length
        dc = float((dens[v] - dens[u]) / (abs(dens[u]) + 1e-6))
        senders.append(u)
        receivers.append(v)
        e_feat.append([length, float(ddir[0]), float(ddir[1]), float(ddir[2]), dc])
        if not args.no_bidirectional_edges:
            senders.append(v)
            receivers.append(u)
            e_feat.append([length, float(-ddir[0]), float(-ddir[1]), float(-ddir[2]), -dc])

    senders = np.asarray(senders, dtype=np.int32)
    receivers = np.asarray(receivers, dtype=np.int32)
    edges = np.asarray(e_feat, dtype=np.float32)
    print(f"Directed edges: {senders.size:,}, edge_feat_dim={edges.shape[1] if edges.size else 0}")

    graph = jraph.GraphsTuple(
        nodes=nodes,
        edges=edges,
        senders=senders,
        receivers=receivers,
        n_node=np.asarray([N], dtype=np.int32),
        n_edge=np.asarray([senders.size], dtype=np.int32),
        globals=None,
    )

    payload = {
        "graph": graph,
        "node_feature_names": list(NODE_FEATURE_NAMES),
        "galaxy_table": gal_sub,
        "xyz_mpc": xyz,
        "coordinate_units": "Mpc",
        "wedge_bounds": {
            "ra_min": float(args.ra_min),
            "ra_max": float(args.ra_max),
            "dec_min": float(args.dec_min),
            "dec_max": float(args.dec_max),
            "z_min": float(args.z_min),
            "z_max": float(args.z_max),
            "z_col": str(args.use_z_col),
        },
        "source": {
            "graph_path": str(graph_path),
            "features_path": str(features_path),
            "zcat_path": str(zcat_path),
            "repo_root": str(repo),
        },
    }

    out = args.out_cache.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Wrote DESI wedge Jraph cache: {out}")


if __name__ == "__main__":
    main()

