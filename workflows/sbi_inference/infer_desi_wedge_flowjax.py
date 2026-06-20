#!/usr/bin/env python3
"""Run an Abacus-trained FlowJAX NPE (GNN encoder + normalizing flow) on a real
DESI LOA wedge, producing PER-GALAXY posterior eigenvalue samples and T-Web class
probabilities. This is the SBI/NPE companion to
``jraph_infer_desi_wedge_from_gnn_npz.py`` (which runs the 15-d Jraph *regression*).

Unlike the regression path, the NPE yields a full posterior per galaxy, so we get
posterior-mean λ, posterior WIDTH (uncertainty — the regression has none) and
class PROBABILITIES (not just a hard class). DESI has no per-galaxy truth, so the
truth-requiring calibration (TARP/SBC) lives on the Abacus side; here we run the
truth-free checks (class-prob count/marginal consistency, sensible class fractions).

Parity (must match what the GNN encoder saw in training — verified against
``build_abacus_sbi_cache._build_graph_from_npz``):
  * nodes: ``node_feature_scaler.transform(x_raw + 1e-6)`` (box-cox from the cache)
  * edges: bidirectional dup + log(edge_length,density_contrast) + StandardScaler,
    where the scaler is fit on the **path1 fiberassign** Abacus wedge npz (the
    FlowJAX model's actual training graph) — NOT the regression wedge. Both wedges
    have 100935 nodes, so a shape check can't catch a mix-up; we hard-assert the
    fitted scaler constants instead.

Outputs (in --output-dir/<run_name>/):
  - desi_wedge_flowjax_preds.npz  (global_node_id, ra, dec, z, lambda_mean[N,3],
                                   lambda_std[N,3], classprob[N,4], hard_class[N],
                                   p_exceed[N,3])
  - summary.json                  (provenance, parity constants, class fractions,
                                   reference fractions for the comparison plot)
Figures are produced by the companion plot_desi_wedge_flowjax.py (no GPU needed).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
# NOTE: do NOT force JAX_PLATFORMS=cpu — this run uses the GPU.

# shared.* must resolve to the Illustris training repo, which also defines the
# FlowJAX model/encoder. GraphWeb_DESI also has a top-level shared/ — put Illustris
# FIRST so `import shared.*` and `from workflows.sbi...` resolve there.
_ILLUSTRIS_ROOT = Path(
    os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")
).expanduser().resolve()
if not (_ILLUSTRIS_ROOT / "shared" / "graph_net_models.py").exists():
    raise FileNotFoundError(f"Illustris repo not found at {_ILLUSTRIS_ROOT}; set ILLUSTRIS_ROOT.")
if str(_ILLUSTRIS_ROOT) not in sys.path:
    sys.path.insert(0, str(_ILLUSTRIS_ROOT))

import jax
import jraph  # noqa: F401  (kept for parity with the regression script's deps)

from workflows.sbi.plot_flowjax_posteriors import (
    load_flowjax_model,
    create_gnn_and_flow,
    batched_sample_posterior,
)
from shared.eigenvalue_transformations import samples_to_raw_eigenvalues, posterior_to_classprobs

_GRAPHWEB_ROOT = Path(__file__).resolve().parents[2]

NODE_FEATURE_NAMES = ("Degree", "Clustering", "Density", "Neigh Density", "I_eig1", "I_eig2", "I_eig3")
CLASS_ORDER = ("void", "wall", "filament", "cluster")

# Edge-scaler constants for the path1 fiberassign wedge (verified). If the wrong
# --abacus-gnn-arrays is passed (e.g. the regression wedge), these differ -> abort.
_PATH1_EDGE_SCALER_MEAN = np.array([2.1516, 0.0])
_PATH1_EDGE_SCALER_SCALE = np.array([0.6963, 1.7736])


def _load_abacus_gnn_parity():
    path = _GRAPHWEB_ROOT / "shared" / "abacus_gnn_parity.py"
    spec = importlib.util.spec_from_file_location("abacus_gnn_parity", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main(args):
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = (args.output_dir / args.run_name).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    parity = _load_abacus_gnn_parity()

    print(f"JAX devices: {jax.devices()}", flush=True)

    # --- model + flow ---
    model_path = args.model_path.expanduser().resolve()
    gnn_params, config, target_scaler, flow_filename, increment_mode = load_flowjax_model(str(model_path))
    print(f"increment_mode={increment_mode}; config latent={config['latent_size']} "
          f"passes={config['num_passes']} heads={config['num_heads']}", flush=True)

    # --- calibration cache: node_feature_scaler (+ Abacus truth eigenvalues for ref) ---
    with args.calibration_cache.expanduser().resolve().open("rb") as f:
        cal = pickle.load(f)
    node_feature_scaler = cal.get("node_feature_scaler")
    if node_feature_scaler is None:
        raise KeyError("Calibration cache has no node_feature_scaler; the FlowJAX cache was built "
                       "with --power-scale-node-features, so this must be present.")
    abacus_true_eigs = (np.asarray(cal["eigenvalues_raw"], dtype=np.float64)
                        if cal.get("eigenvalues_raw") is not None else None)

    # --- DESI wedge metadata + coordinate-units guard ---
    with args.desi_gnn_metadata.expanduser().resolve().open() as f:
        gmeta = json.load(f)
    coord_units = str(gmeta.get("coordinate_units", "unknown"))
    if coord_units.lower() in ("mpc/h", "mpc_per_h", "unknown"):
        raise SystemExit(f"DESI coordinate_units={coord_units!r}; Abacus parity requires 'mpc'. Rebuild wedge.")
    node_cols = list(gmeta.get("node_feature_columns", []))
    if node_cols and tuple(node_cols) != NODE_FEATURE_NAMES:
        raise ValueError(f"DESI node_feature_columns mismatch:\n exp {NODE_FEATURE_NAMES}\n got {tuple(node_cols)}")

    # --- DESI wedge arrays ---
    desi = np.load(args.desi_gnn_arrays.expanduser().resolve())
    x_raw = np.asarray(desi["x"], dtype=np.float64)
    edge_index = np.asarray(desi["edge_index"], dtype=np.int64)
    edge_attr = np.asarray(desi["edge_attr"], dtype=np.float32)
    if x_raw.ndim != 2 or x_raw.shape[1] != 7:
        raise ValueError(f"DESI x must be (N,7); got {x_raw.shape}")

    # --- edges: fit scaler on PATH1 wedge, assert constants, transform DESI edges ---
    edge_scaler = parity.fit_edge_length_density_scaler_from_gnn_npz(
        args.abacus_gnn_arrays.expanduser().resolve(), make_bidirectional=True)
    if not (np.allclose(edge_scaler.mean_, _PATH1_EDGE_SCALER_MEAN, atol=0.05)
            and np.allclose(edge_scaler.scale_, _PATH1_EDGE_SCALER_SCALE, atol=0.05)):
        raise SystemExit(
            "Edge-scaler constants do not match the path1 fiberassign wedge — wrong --abacus-gnn-arrays?\n"
            f"  got   mean={edge_scaler.mean_}, scale={edge_scaler.scale_}\n"
            f"  want  mean≈{_PATH1_EDGE_SCALER_MEAN}, scale≈{_PATH1_EDGE_SCALER_SCALE}")
    edge_index_b, edge_attr_b = parity.prepare_edges_for_jraph_forward(
        edge_index, edge_attr, edge_scaler, make_bidirectional=True)
    print(f"edge scaler OK: mean={edge_scaler.mean_.round(4)} scale={edge_scaler.scale_.round(4)}; "
          f"bidir edges={edge_attr_b.shape[0]}", flush=True)

    # Phase-0a domain correction: the Abacus-fit scaler leaves DESI's scaled edge
    # features off the training N(0,1) (graph-scale shift, ~-0.11σ on edge_length).
    # Re-standardise the scaled edge_length (col 0) and density_contrast (col 4) to
    # zero-mean/unit-std so the GNN sees the training edge distribution. Tests whether
    # the cluster deficit is caused by the edge-scale domain shift.
    if args.edge_domain_adapt:
        edge_attr_b = np.array(edge_attr_b, dtype=np.float32, copy=True)
        for col in (0, 4):
            m, s = float(edge_attr_b[:, col].mean()), float(edge_attr_b[:, col].std())
            edge_attr_b[:, col] = (edge_attr_b[:, col] - m) / (s + 1e-9)
            print(f"[edge-domain-adapt] col{col}: DESI mean {m:+.3f} std {s:.3f} -> N(0,1)", flush=True)

    # --- nodes: box-cox (match training) + distribution parity check ---
    x = node_feature_scaler.transform(x_raw + 1e-6).astype(np.float32)
    col_mean, col_std = x.mean(0), x.std(0)
    print("node post-box-cox per-col mean:", col_mean.round(3), flush=True)
    print("node post-box-cox per-col std :", col_std.round(3), flush=True)
    if (np.abs(col_mean) > 0.5).any() or ((col_std < 0.7) | (col_std > 1.4)).any():
        print("WARNING: DESI node distribution drifted from training (~0/1); transfer may be degraded.", flush=True)

    n_nodes, n_edges = int(x.shape[0]), int(edge_index_b.shape[1])
    graph = jraph.GraphsTuple(
        nodes=x, edges=edge_attr_b,
        senders=edge_index_b[0].astype(np.int32), receivers=edge_index_b[1].astype(np.int32),
        n_node=np.asarray([n_nodes], dtype=np.int32), n_edge=np.asarray([n_edges], dtype=np.int32),
        globals=None)

    # --- forward: GNN encoder -> embeddings -> posterior samples -> raw eigenvalues ---
    gnn, flow = create_gnn_and_flow(config, flow_filename, graph, jax.random.key(args.seed))
    embeddings = np.asarray(gnn.apply(gnn_params, jax.random.key(0), graph, is_training=False))
    if embeddings.shape != (n_nodes, config["latent_size"]):
        raise ValueError(f"embeddings shape {embeddings.shape} != ({n_nodes},{config['latent_size']})")
    print(f"embeddings: {embeddings.shape}; sampling {args.num_posterior_samples}/galaxy ...", flush=True)

    samples_scaled = batched_sample_posterior(
        flow, embeddings, args.num_posterior_samples, jax.random.key(123), chunk_size=args.chunk_size)
    samples_raw = samples_to_raw_eigenvalues(samples_scaled, target_scaler, increment_mode)  # [N,K,3]

    # TRUE ordering-violation rate (no ground truth needed): fraction of the flow's
    # own posterior samples whose inverted (λ1,λ2,λ3) are NOT ascending — i.e. the
    # flow drew a negative linear increment. This is a self-consistency property of
    # the model's samples, not an accuracy claim.
    asc = (samples_raw[..., 0] <= samples_raw[..., 1]) & (samples_raw[..., 1] <= samples_raw[..., 2])
    viol_rate = float(np.mean(~asc))
    print(f"TRUE ordering-violation rate (non-ascending samples): {viol_rate:.4f}", flush=True)
    # Post-hoc fix: sort each posterior sample ascending. Physical eigenvalues are
    # ordered and the T-web count is order-independent, so this is a monotonic,
    # no-retraining correction applied at eval time.
    if not args.no_sort:
        samples_raw = np.sort(samples_raw, axis=-1)

    # --- per-galaxy posterior products ---
    lambda_mean = samples_raw.mean(axis=1)                # [N,3]
    lambda_std = samples_raw.std(axis=1)                  # [N,3] (NPE uncertainty)
    cp = posterior_to_classprobs(samples_raw, lambda_th=args.lambda_threshold)
    classprob = np.stack([np.asarray(cp[c]) for c in CLASS_ORDER], axis=1)  # [N,4]
    hard_class = classprob.argmax(axis=1).astype(np.int8)
    consistency = float(cp["consistency_max_abs_diff"])
    print(f"class-prob count/marginal consistency: {consistency:.3e}", flush=True)
    if consistency >= 1e-3:
        print("WARNING: high consistency value -> flow emitted out-of-order λ samples.", flush=True)

    npe_fracs = np.array([float(np.mean(hard_class == k)) for k in range(4)])
    print("DESI NPE class fractions (argmax):", {CLASS_ORDER[k]: round(npe_fracs[k], 3) for k in range(4)}, flush=True)
    if npe_fracs.max() > 0.85:
        print("WARNING: a class exceeds 0.85 — possible scaler mismatch / collapse.", flush=True)

    # --- ra/dec/z + global ids for the sky maps and downstream joins ---
    global_ids = np.load(args.desi_global_node_ids.expanduser().resolve()).astype(np.int64)
    wc = np.load(args.desi_wedge_catalog_npz.expanduser().resolve())
    ra, dec, zz = (np.asarray(wc["ra"], np.float64), np.asarray(wc["dec"], np.float64), np.asarray(wc["z"], np.float64))

    # Abacus truth class fractions at this threshold (for the comparison plot).
    abacus_fracs = None
    if abacus_true_eigs is not None:
        acp = posterior_to_classprobs(abacus_true_eigs[:, None, :], lambda_th=args.lambda_threshold)
        abacus_fracs = {c: float(np.mean(acp[c])) for c in CLASS_ORDER}

    # Subsample of full posterior samples (for KDE distributions / the skewer) — full
    # [N,128,3] would be ~170 MB; keep a representative random subset instead.
    rng_sub = np.random.default_rng(0)
    sub = np.sort(rng_sub.choice(n_nodes, min(n_nodes, args.save_sample_subset), replace=False))
    np.savez_compressed(
        out_dir / "desi_wedge_flowjax_preds.npz",
        global_node_id=global_ids, ra=ra, dec=dec, z=zz,
        lambda_mean=lambda_mean, lambda_std=lambda_std,
        classprob=classprob, hard_class=hard_class, p_exceed=np.asarray(cp["p_exceed"]),
        embeddings=embeddings.astype(np.float32),
        sample_subset_idx=sub, lambda_samples_subset=samples_raw[sub].astype(np.float32))
    print(f"Saved: {out_dir/'desi_wedge_flowjax_preds.npz'}", flush=True)

    summary = {
        "timestamp": ts,
        "model_path": str(model_path),
        "flow_filename": str(flow_filename),
        "calibration_cache": str(args.calibration_cache),
        "abacus_gnn_arrays_for_edge_scaler": str(args.abacus_gnn_arrays),
        "increment_mode": increment_mode,
        "lambda_threshold": float(args.lambda_threshold),
        "num_posterior_samples": int(args.num_posterior_samples),
        "node_feature_transform": "box-cox PowerTransformer(+1e-6) from calibration cache",
        "edge_attr_transform": "bidirectional + log(len,contrast) + StandardScaler (path1 wedge fit)",
        "edge_scaler_mean": edge_scaler.mean_.tolist(),
        "edge_scaler_scale": edge_scaler.scale_.tolist(),
        "node_post_boxcox_mean": col_mean.tolist(),
        "node_post_boxcox_std": col_std.tolist(),
        "classprob_consistency_max_abs_diff": consistency,
        "true_ordering_violation_rate": viol_rate,
        "post_hoc_sort_applied": (not args.no_sort),
        "desi": {"n_nodes": n_nodes, "n_edges": n_edges,
                 "class_fractions_npe": {CLASS_ORDER[k]: float(npe_fracs[k]) for k in range(4)}},
        "reference_fractions": {
            "abacus_truth": abacus_fracs,
            "regression_desi": {"void": 0.246, "wall": 0.462, "filament": 0.266, "cluster": 0.026},
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"Saved: {out_dir/'summary.json'}\nDone.", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="FlowJAX NPE inference on a DESI wedge")
    ap.add_argument("--model-path", type=Path, required=True)
    ap.add_argument("--calibration-cache", type=Path, required=True,
                    help="FlowJAX training cache (node_feature_scaler + eigenvalues_raw)")
    ap.add_argument("--abacus-gnn-arrays", type=Path, required=True,
                    help="PATH1 fiberassign wedge raw cugraph npz (edge scaler fit); asserted")
    ap.add_argument("--desi-gnn-arrays", type=Path, required=True)
    ap.add_argument("--desi-gnn-metadata", type=Path, required=True)
    ap.add_argument("--desi-global-node-ids", type=Path, required=True)
    ap.add_argument("--desi-wedge-catalog-npz", type=Path, required=True)
    ap.add_argument("--num-posterior-samples", type=int, default=128)
    ap.add_argument("--lambda-threshold", type=float, default=0.2)
    ap.add_argument("--chunk-size", type=int, default=512)
    ap.add_argument("--edge-domain-adapt", action="store_true",
                    help="Phase-0a: re-standardise DESI scaled edge_length+density_contrast to "
                         "training N(0,1) (corrects the graph-scale domain shift).")
    ap.add_argument("--no-sort", action="store_true",
                    help="Disable the post-hoc ascending sort of posterior samples (keep raw flow order).")
    ap.add_argument("--save-sample-subset", type=int, default=20000,
                    help="How many galaxies' full posterior samples to save (for KDE/skewer).")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output-dir", type=Path, default=Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs"))
    ap.add_argument("--run-name", type=str, default="desi_wedge_flowjax_linear")
    main(ap.parse_args())
