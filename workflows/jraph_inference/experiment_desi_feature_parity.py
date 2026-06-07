#!/usr/bin/env python3
"""One-off parity experiments: Mpc/h scaling, node z-score, Abacus edge scaling."""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import jax
import haiku as hk
import jraph
import numpy as np
from sklearn.preprocessing import StandardScaler

ILL = Path("/global/homes/d/dkololgi/TNG/Illustris")
sys.path.insert(0, str(ILL))
from shared.graph_net_models import make_graph_network  # noqa: E402

CACHE = Path(
    "/pscratch/sd/d/dkololgi/abacus/sbi_caches/"
    "staged_mock_stage3_postcollision_full_rs7_wedge_v_limited_incomplete_expanded_"
    "ra120_160_dec14p5_30p6_z0p2_0p3_sbi_cache_15d_halo_xcom.pkl"
)
ABACUS_GNN = Path(
    "/pscratch/sd/d/dkololgi/abacus/graph_constructions/"
    "staged_mock_stage3_postcollision_full_rs7_wedge_v_limited_incomplete_expanded_"
    "ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz"
)
DESI_GNN = Path(
    "/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_from_fullgraph/"
    "desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_gnn_arrays.npz"
)
DESI_XYZ = Path(
    "/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_from_fullgraph/"
    "desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_points_xyz.npy"
)
CKPT_JSON = Path(
    "/pscratch/sd/d/dkololgi/abacus/jraph_runs/"
    "staged_mock_stage3_postcollision_full_rs7_wedge_v_limited_incomplete_expanded_"
    "ra120_160_dec14p5_30p6_z0p2_0p3_15d_halo_xcom/checkpoints/best_checkpoint.json"
)
OUT = Path(
    "/pscratch/sd/d/dkololgi/graphweb_desi/inference_outputs/experiment_feature_parity_20260529"
)
H = 0.6766  # Planck18


def eig_from_15(raw: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    l1 = raw[:, 0]
    l2 = l1 + np.maximum(raw[:, 1], eps)
    l3 = l2 + np.maximum(raw[:, 2], eps)
    return np.stack([l1, l2, l3], axis=-1)


def cweb_fracs(lam: np.ndarray, thr: float = 0.2) -> dict[str, float]:
    cls = (lam > thr).sum(axis=1)
    b = np.bincount(cls.astype(np.int64), minlength=4).astype(np.float64)
    b /= max(1, cls.size)
    return {k: float(b[i]) for i, k in enumerate(["void", "wall", "filament", "cluster"])}


def scale_edges_like_abacus_cache(edge_attr: np.ndarray, scaler: StandardScaler) -> np.ndarray:
    e = np.asarray(edge_attr, dtype=np.float32).copy()
    e[:, 0] = np.log(np.maximum(e[:, 0], 1e-6))
    e[:, 4] = np.log(np.maximum(e[:, 4], 1e-6))
    e[:, [0, 4]] = scaler.transform(e[:, [0, 4]]).astype(np.float32)
    return e


def fit_abacus_edge_scaler() -> StandardScaler:
    """Match build_abacus_sbi_cache: bidirectional + log + StandardScaler on cols 0,4."""
    d = np.load(ABACUS_GNN)
    edge_attr = np.asarray(d["edge_attr"], dtype=np.float64)
    edge_index = np.asarray(d["edge_index"], dtype=np.int64)
    rev = edge_attr.copy()
    rev[:, 1:4] *= -1.0
    rev[:, 4] = 1.0 / np.maximum(rev[:, 4], 1e-6)
    edge_attr = np.concatenate([edge_attr, rev], axis=0)
    edge_attr[:, 0] = np.log(np.maximum(edge_attr[:, 0], 1e-6))
    edge_attr[:, 4] = np.log(np.maximum(edge_attr[:, 4], 1e-6))
    sc = StandardScaler()
    sc.fit(edge_attr[:, [0, 4]])
    return sc


def make_graph(x: np.ndarray, edge_index: np.ndarray, edge_attr: np.ndarray) -> jraph.GraphsTuple:
    return jraph.GraphsTuple(
        nodes=np.asarray(x, dtype=np.float32),
        edges=np.asarray(edge_attr, dtype=np.float32),
        senders=edge_index[0].astype(np.int32),
        receivers=edge_index[1].astype(np.int32),
        n_node=np.asarray([x.shape[0]], dtype=np.int32),
        n_edge=np.asarray([edge_index.shape[1]], dtype=np.int32),
        globals=None,
    )


def zscore_nodes_to_abacus(x: np.ndarray, abacus_nodes: np.ndarray) -> np.ndarray:
    mu = abacus_nodes.mean(axis=0)
    sd = abacus_nodes.std(axis=0)
    sd = np.where(sd < 1e-12, 1.0, sd)
    return ((x - mu) / sd).astype(np.float32)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    from astropy.cosmology import Planck18 as cosmo

    h = float(cosmo.h)
    report: dict = {"h": h}

    # --- Task 1: Mpc/h ---
    ab_raw = np.load(ABACUS_GNN)
    desi_raw = np.load(DESI_GNN)
    ab_el_raw = np.asarray(ab_raw["edge_attr"][:, 0], dtype=np.float64)
    desi_el_raw = np.asarray(desi_raw["edge_attr"][:, 0], dtype=np.float64)
    xyz = np.load(DESI_XYZ)
    ei = np.asarray(desi_raw["edge_index"], dtype=np.int64)
    rng = np.random.default_rng(0)
    sub = rng.choice(ei.shape[1], size=min(100_000, ei.shape[1]), replace=False)
    el_re = np.linalg.norm(xyz[ei[0, sub]] - xyz[ei[1, sub]], axis=1)
    el_mpc = np.linalg.norm(xyz[ei[0, sub]] / h - xyz[ei[1, sub]] / h, axis=1)

    report["task1_mpc_units"] = {
        "abacus_raw_edge_length_mean_mpc": float(ab_el_raw.mean()),
        "abacus_raw_edge_length_median_mpc": float(np.median(ab_el_raw)),
        "desi_raw_edge_length_mean_mpc_per_h": float(desi_el_raw.mean()),
        "desi_recomputed_from_xyz_mean_mpc_per_h": float(el_re.mean()),
        "desi_coords_div_h_mean_mpc": float(el_mpc.mean()),
        "ratio_desi_mpc_per_h_over_abacus_mpc": float(desi_el_raw.mean() / max(ab_el_raw.mean(), 1e-12)),
        "ratio_after_coord_div_h_over_abacus": float(el_mpc.mean() / max(ab_el_raw.mean(), 1e-12)),
        "note": "Abacus export uses comoving Mpc; DESI builder uses Mpc/h (×h). Pure unit fix gives edge ratio ~1/h, not ~70M.",
    }

    with open(CACHE, "rb") as f:
        cal = pickle.load(f)
    ab_nodes = np.asarray(cal["graph"].nodes, dtype=np.float64)
    scaler_t = cal["target_scaler"]
    edge_scaler = fit_abacus_edge_scaler()

    desi_x = np.asarray(desi_raw["x"], dtype=np.float32)
    desi_ea = np.asarray(desi_raw["edge_attr"], dtype=np.float32)
    desi_ei = np.asarray(desi_raw["edge_index"], dtype=np.int64)

    with open(CKPT_JSON) as f:
        model_path = Path(json.load(f)["path"])
    with open(model_path, "rb") as f:
        params = pickle.load(f)["params"]

    net = hk.transform(
        make_graph_network(num_passes=8, latent_size=96, num_heads=8, dropout_rate=0.15, output_dim=15)
    )

    @jax.jit
    def predict(p, g, rng):
        return net.apply(p, rng, g, is_training=False).nodes

    variants = {
        "baseline_unmodified": (desi_x, desi_ea),
        "nodes_zscore_abacus_moments": (zscore_nodes_to_abacus(desi_x, ab_nodes), desi_ea),
        "edges_abacus_log_standardize": (desi_x, scale_edges_like_abacus_cache(desi_ea, edge_scaler)),
        "nodes_zscore_and_edges_scaled": (
            zscore_nodes_to_abacus(desi_x, ab_nodes),
            scale_edges_like_abacus_cache(desi_ea, edge_scaler),
        ),
    }

    report["task2_forward"] = {}
    rng = jax.random.PRNGKey(999)
    for name, (xn, ea) in variants.items():
        g = make_graph(xn, desi_ei, ea)
        raw15 = scaler_t.inverse_transform(np.asarray(predict(params, g, rng)))
        lam = eig_from_15(raw15)
        report["task2_forward"][name] = {
            "class_fractions_thr0p2": cweb_fracs(lam),
            "raw_ch0_2_mean": raw15[:, :3].mean(axis=0).tolist(),
            "frac_delta2_lt_0": float((raw15[:, 1] < 0).mean()),
            "frac_delta3_lt_0": float((raw15[:, 2] < 0).mean()),
        }

    # Baseline from saved eigfix
    eigfix = np.load(
        "/pscratch/sd/d/dkololgi/graphweb_desi/inference_outputs/"
        "infer_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_best_epoch25427_eigfix/preds_raw_15d.npy"
    )
    raw_base = scaler_t.inverse_transform(
        np.asarray(predict(params, make_graph(desi_x, desi_ei, desi_ea), jax.random.PRNGKey(999)))
    )
    report["task2_forward"]["saved_eigfix_baseline"] = {
        "class_fractions_thr0p2": cweb_fracs(eig_from_15(eigfix)),
        "max_abs_diff_raw15_vs_rerun": float(np.max(np.abs(eigfix - raw_base))),
    }

    out_json = OUT / "experiment_feature_parity_report.json"
    out_json.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    print("Wrote", out_json)


if __name__ == "__main__":
    main()
