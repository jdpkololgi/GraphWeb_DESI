#!/usr/bin/env python3
"""Run an Abacus-trained Jraph regression checkpoint on a DESI wedge (GNN arrays).

This is the inference companion to the Abacus regression training/eval pipeline:
  - model checkpoint: from an Abacus jraph run directory (e.g. .../checkpoints/*best.pkl)
  - calibration cache: Abacus wedge SBI cache used for training (provides target_scaler and true eigenvalues)
  - DESI wedge inputs: output of `subset_desi_graph_wedge.py` (wedge-specific `*_gnn_arrays.npz`
    and small wedge catalog arrays for diagnostics)

Outputs (in --output-dir/<run_name>/):
  - preds_scaled_15d.npy
  - preds_raw_15d.npy                 (inverse-transformed by Abacus target_scaler)
  - preds_lambda123.npy               (physical λ1,λ2,λ3 reconstructed from increments; for plots/classes)
  - class_lambda_thr0p2.npy           (0..3 void/wall/filament/cluster by lambda_threshold)
  - summary.json                      (paths + class fractions)
  - plots: eigenvalue histograms and class-fraction bar chart

Note on eigenvalues (15-d halo_xcom caches):
After ``target_scaler.inverse_transform``, the first three channels are **ordered linear
increments** (same as ``build_abacus_sbi_cache.py`` with derivative columns):
  v1 = λ1,  v2 = λ2 − λ1,  v3 = λ3 − λ2.
Physical eigenvalues are reconstructed as λ2 = λ1 + max(v2, ε), λ3 = λ2 + max(v3, ε).
Do **not** treat v2/v3 as λ2/λ3 directly. (Softplus increments apply only to 3-d caches without
derivatives and ``use_transformed_eig=True`` — not this 15-d run.)

Node features at forward: if the calibration cache contains ``node_feature_scaler``
(box-cox ``PowerTransformer`` from ``build_abacus_sbi_cache.py --power-scale-node-features``),
apply ``scaler.transform(x + 1e-6)`` to DESI raw ``x`` before the GNN forward pass.
Training/eval loads the **already transformed** graph from the SBI cache; feeding raw DESI
features without this step collapses predictions (e.g. ~100% void).

Edge attributes at forward: duplicate bidirectional edges (Abacus cache convention), then
``log(edge_length)`` / ``log(density_contrast)`` and ``StandardScaler`` fit on the Abacus wedge
GNN NPZ. Pass ``--abacus-gnn-arrays``. GNN NPZ files store physical edge features only.
"""

from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

# JAX: avoid CUDA plugin initialization on CPU-only nodes (otherwise jax prints
# CUDA_ERROR_NO_DEVICE / xla_cuda12.initialize errors even though CPU fallback works).
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
os.environ.setdefault("PYTHONNOUSERSITE", "1")
if not os.environ.get("JAX_PLATFORMS", "").strip():
    os.environ["JAX_PLATFORMS"] = "cpu"
if os.environ.get("JAX_PLATFORMS", "cpu").strip().lower() == "cpu":
    # Prevent PJRT from probing CUDA when we only use CPU (typical on Perlmutter CPU interactive nodes).
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

# `shared.graph_net_models` lives in the Illustris training repo (same as jraph_pipeline),
# not under GraphWeb_DESI.
# GraphWeb_DESI also has a top-level `shared/` package — do NOT put GraphWeb_DESI on sys.path
# before Illustris, or `import shared.*` resolves to the wrong package.
_ILLUSTRIS_ROOT = Path(
    os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")
).expanduser().resolve()
if _ILLUSTRIS_ROOT.is_dir() and str(_ILLUSTRIS_ROOT) not in sys.path:
    sys.path.insert(0, str(_ILLUSTRIS_ROOT))
if not (_ILLUSTRIS_ROOT / "shared" / "graph_net_models.py").exists():
    raise FileNotFoundError(
        f"Illustris repo not found at {_ILLUSTRIS_ROOT}. Set ILLUSTRIS_ROOT to the directory "
        "that contains shared/graph_net_models.py (same tree as workflows/jraph/)."
    )

import importlib.util

import jax
import haiku as hk
import jraph
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from shared.graph_net_models import make_graph_network

_GRAPHWEB_ROOT = Path(__file__).resolve().parents[2]


def _load_abacus_gnn_parity():
    path = _GRAPHWEB_ROOT / "shared" / "abacus_gnn_parity.py"
    spec = importlib.util.spec_from_file_location("abacus_gnn_parity", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_abacus_parity = _load_abacus_gnn_parity()


NODE_FEATURE_NAMES = (
    "Degree",
    "Clustering",
    "Density",
    "Neigh Density",
    "I_eig1",
    "I_eig2",
    "I_eig3",
)


def _load_params(model_path: Path) -> tuple[object, int | None]:
    with model_path.open("rb") as f:
        obj = pickle.load(f)
    if isinstance(obj, dict) and "params" in obj:
        return obj["params"], obj.get("epoch")
    return obj, None


def _eig_from_15d_linear(raw15: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    """Reconstruct (λ1, λ2, λ3) from 15-d linear increment targets (channels 0–2)."""
    raw15 = np.asarray(raw15, dtype=np.float64)
    l1 = raw15[:, 0]
    d2 = np.maximum(raw15[:, 1], eps)
    d3 = np.maximum(raw15[:, 2], eps)
    l2 = l1 + d2
    l3 = l2 + d3
    return np.stack([l1, l2, l3], axis=-1)


def _lambda_threshold_classes(lam123: np.ndarray, thr: float) -> np.ndarray:
    """Return T-web-like class by number of eigenvalues > thr (0..3)."""
    lam = np.asarray(lam123, dtype=np.float64)
    return (lam > float(thr)).sum(axis=1).astype(np.int8)


def _hist_compare(a: np.ndarray, b: np.ndarray, out_png: Path, title: str, labels: tuple[str, str]) -> None:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.0))
    bins = 80
    ax.hist(a, bins=bins, density=True, alpha=0.55, label=labels[0])
    ax.hist(b, bins=bins, density=True, alpha=0.55, label=labels[1])
    ax.set_xlabel("eigenvalue")
    ax.set_ylabel("density")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _bar_class_fractions(fracs: np.ndarray, out_png: Path, title: str) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(6.0, 3.8))
    names = ["void (0)", "wall (1)", "filament (2)", "cluster (3)"]
    ax.bar(np.arange(4), fracs, color="tab:blue", alpha=0.85)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.set_ylabel("fraction")
    ax.set_ylim(0.0, 1.0)
    ax.set_title(title)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)

    p.add_argument(
        "--abacus-run-dir",
        type=Path,
        required=True,
        help="Abacus wedge regression run directory (contains checkpoints/best_checkpoint.json).",
    )
    p.add_argument(
        "--calibration-cache",
        type=Path,
        required=True,
        help="Abacus wedge SBI cache used for training (provides target_scaler and true eigenvalues).",
    )
    p.add_argument(
        "--abacus-gnn-arrays",
        type=Path,
        required=True,
        help="Abacus wedge *_cugraph_gnn_arrays.npz (fits edge log+StandardScaler; must match training wedge).",
    )
    p.add_argument(
        "--desi-gnn-arrays",
        type=Path,
        required=True,
        help="DESI wedge *_gnn_arrays.npz (from subset_desi_graph_wedge.py).",
    )
    p.add_argument(
        "--desi-gnn-metadata",
        type=Path,
        required=True,
        help="DESI wedge *_gnn_metadata.json (from subset_desi_graph_wedge.py).",
    )
    p.add_argument(
        "--desi-global-node-ids",
        type=Path,
        required=True,
        help="DESI wedge *_global_node_ids.npy (maps local wedge -> global catalog row id).",
    )
    p.add_argument(
        "--desi-wedge-catalog-npz",
        type=Path,
        required=True,
        help="DESI wedge *_wedge_catalog_minimal.npz (ra/dec/z arrays aligned to wedge node order).",
    )

    p.add_argument("--output-dir", type=Path, required=True, help="Parent directory for this inference run.")
    p.add_argument("--run-name", type=str, default=None, help="Optional run folder name. Default: infer_<timestamp>.")

    # Model hparams (default to your wedge run values).
    p.add_argument("--latent-size", type=int, default=96)
    p.add_argument("--num-heads", type=int, default=8)
    p.add_argument("--num-passes", type=int, default=8)
    p.add_argument("--dropout", type=float, default=0.15)
    p.add_argument("--seed", type=int, default=42)

    p.add_argument("--lambda-threshold", type=float, default=0.2)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    abacus_run = args.abacus_run_dir.expanduser().resolve()
    best_json = abacus_run / "checkpoints" / "best_checkpoint.json"
    if not best_json.exists():
        raise FileNotFoundError(f"Missing {best_json}; expected an Abacus jraph run dir.")
    best = json.loads(best_json.read_text(encoding="utf-8"))
    model_path = Path(best["path"]).expanduser().resolve()
    if not model_path.exists():
        raise FileNotFoundError(f"Best checkpoint path not found: {model_path}")

    out_root = args.output_dir.expanduser().resolve()
    run_name = args.run_name or f"infer_{ts}"
    out_dir = out_root / run_name
    out_dir.mkdir(parents=True, exist_ok=True)

    gmeta = json.loads(args.desi_gnn_metadata.expanduser().resolve().read_text(encoding="utf-8"))
    coord_units = gmeta.get("coordinate_units", "unknown")
    if str(coord_units).lower() in ("mpc/h", "mpc_per_h", "unknown"):
        print(
            f"WARNING: DESI gnn metadata coordinate_units={coord_units!r}; "
            "Abacus training uses Mpc. Rebuild graph with --coord-units mpc.",
            flush=True,
        )

    # Load calibration cache early (node scaler + target_scaler).
    cal_path = args.calibration_cache.expanduser().resolve()
    with cal_path.open("rb") as f:
        cal = pickle.load(f)
    target_scaler = cal.get("target_scaler")
    if target_scaler is None:
        raise KeyError("Calibration cache missing target_scaler.")
    node_feature_scaler = cal.get("node_feature_scaler")
    abacus_true_eigs = (
        np.asarray(cal.get("eigenvalues_raw"), dtype=np.float64)
        if cal.get("eigenvalues_raw") is not None
        else None
    )

    # Load DESI wedge GNN arrays (physical / unscaled node features in NPZ).
    desi_npz = np.load(args.desi_gnn_arrays.expanduser().resolve())
    x_raw = np.asarray(desi_npz["x"], dtype=np.float64)
    edge_index = np.asarray(desi_npz["edge_index"], dtype=np.int64)
    edge_attr = np.asarray(desi_npz["edge_attr"], dtype=np.float32)
    if edge_index.shape[0] != 2:
        raise ValueError(f"edge_index must be (2,E); got {edge_index.shape}")

    edge_scaler = _abacus_parity.fit_edge_length_density_scaler_from_gnn_npz(
        args.abacus_gnn_arrays.expanduser().resolve(),
        make_bidirectional=True,
    )
    edge_index, edge_attr = _abacus_parity.prepare_edges_for_jraph_forward(
        edge_index,
        edge_attr,
        edge_scaler,
        make_bidirectional=True,
    )

    # Sanity-check node feature order matches the Abacus-trained model expectation (7 features).
    node_cols = list(gmeta.get("node_feature_columns", []))
    if node_cols and tuple(node_cols) != NODE_FEATURE_NAMES:
        raise ValueError(f"DESI node_feature_columns mismatch.\nExpected: {NODE_FEATURE_NAMES}\nGot: {tuple(node_cols)}")
    if x_raw.ndim != 2 or x_raw.shape[1] != 7:
        raise ValueError(f"DESI x must be (N,7); got {x_raw.shape}")
    if edge_attr.ndim != 2 or edge_attr.shape[1] != 5:
        raise ValueError(f"DESI edge_attr must be (E,5); got {edge_attr.shape}")

    node_feature_transform = "none"
    if node_feature_scaler is not None:
        x = node_feature_scaler.transform(x_raw + 1e-6).astype(np.float32)
        node_feature_transform = (
            f"box-cox PowerTransformer (+1e-6 shift); method={getattr(node_feature_scaler, 'method', 'box-cox')}"
        )
        print(f"Applied node_feature_scaler from calibration cache ({node_feature_transform})", flush=True)
    else:
        x = x_raw.astype(np.float32)
        print(
            "WARNING: calibration cache has no node_feature_scaler; using raw DESI node features. "
            "This is only valid if training did not use --power-scale-node-features.",
            flush=True,
        )

    n_nodes = int(x.shape[0])
    n_edges = int(edge_index.shape[1])

    # Build a single-graph GraphsTuple.
    senders = edge_index[0].astype(np.int32, copy=False)
    receivers = edge_index[1].astype(np.int32, copy=False)
    graph = jraph.GraphsTuple(
        nodes=x,
        edges=edge_attr,
        senders=senders,
        receivers=receivers,
        n_node=np.asarray([n_nodes], dtype=np.int32),
        n_edge=np.asarray([n_edges], dtype=np.int32),
        globals=None,
    )

    # Load model params and run prediction.
    params, ckpt_epoch = _load_params(model_path)
    net_fn = make_graph_network(
        num_passes=int(args.num_passes),
        latent_size=int(args.latent_size),
        num_heads=int(args.num_heads),
        dropout_rate=float(args.dropout),
        output_dim=15,
    )
    net = hk.transform(net_fn)

    @jax.jit
    def predict(p, g, rng):
        return net.apply(p, rng, g, is_training=False).nodes

    rng = jax.random.PRNGKey(int(args.seed) + 999)
    preds_scaled = np.asarray(predict(params, graph, rng))
    if preds_scaled.ndim != 2 or preds_scaled.shape[1] != 15:
        raise ValueError(f"Unexpected preds shape {preds_scaled.shape}; expected (N,15).")

    preds_raw15 = target_scaler.inverse_transform(preds_scaled)
    preds_raw15 = np.asarray(preds_raw15, dtype=np.float64)
    preds_lam = _eig_from_15d_linear(preds_raw15)

    # Lambda-threshold classes for DESI predictions + Abacus true (if available).
    thr = float(args.lambda_threshold)
    desi_cls = _lambda_threshold_classes(preds_lam, thr)
    desi_fracs = np.bincount(desi_cls.astype(np.int64), minlength=4).astype(np.float64)
    desi_fracs /= max(1.0, float(desi_cls.size))

    # Write arrays.
    np.save(out_dir / "preds_scaled_15d.npy", np.asarray(preds_scaled, dtype=np.float32))
    np.save(out_dir / "preds_raw_15d.npy", preds_raw15)
    np.save(out_dir / "preds_lambda123.npy", preds_lam)
    np.save(out_dir / f"class_lambda_thr{thr:g}.npy", desi_cls)

    # Save minimal provenance + fractions.
    desi_global_ids = np.load(args.desi_global_node_ids.expanduser().resolve()).astype(np.int64, copy=False)
    wedge_cat = np.load(args.desi_wedge_catalog_npz.expanduser().resolve())
    ra = np.asarray(wedge_cat["ra"], dtype=np.float64)
    dec = np.asarray(wedge_cat["dec"], dtype=np.float64)
    zz = np.asarray(wedge_cat["z"], dtype=np.float64)

    summary = {
        "timestamp": ts,
        "abacus_run_dir": str(abacus_run),
        "abacus_best_checkpoint": str(model_path),
        "checkpoint_epoch": ckpt_epoch,
        "calibration_cache": str(cal_path),
        "abacus_gnn_arrays_for_edge_scaler": str(args.abacus_gnn_arrays),
        "desi_coordinate_units": coord_units,
        "node_feature_transform": node_feature_transform,
        "edge_attr_transform": "bidirectional_dup + log(edge_length,density_contrast) + StandardScaler (Abacus wedge fit)",
        "desi_inputs": {
            "gnn_arrays": str(args.desi_gnn_arrays),
            "gnn_metadata": str(args.desi_gnn_metadata),
            "global_node_ids": str(args.desi_global_node_ids),
            "wedge_catalog_npz": str(args.desi_wedge_catalog_npz),
        },
        "hparams": {
            "latent_size": int(args.latent_size),
            "num_heads": int(args.num_heads),
            "num_passes": int(args.num_passes),
            "dropout": float(args.dropout),
            "seed": int(args.seed),
        },
        "lambda_threshold": thr,
        "desi": {
            "n_nodes": int(n_nodes),
            "n_edges": int(n_edges),
            "class_fractions": {
                "void": float(desi_fracs[0]),
                "wall": float(desi_fracs[1]),
                "filament": float(desi_fracs[2]),
                "cluster": float(desi_fracs[3]),
            },
        },
        "outputs": {
            "preds_scaled_15d": str(out_dir / "preds_scaled_15d.npy"),
            "preds_raw_15d": str(out_dir / "preds_raw_15d.npy"),
            "preds_lambda123": str(out_dir / "preds_lambda123.npy"),
        },
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    # Diagnostic plots: DESI predicted λ distributions + Abacus true λ distributions (if available).
    for k, name in enumerate(["lambda1", "lambda2", "lambda3"]):
        out_png = out_dir / f"hist_{name}_desi_pred_vs_abacus_true.png"
        if abacus_true_eigs is not None and abacus_true_eigs.size > 0:
            _hist_compare(
                preds_lam[:, k],
                abacus_true_eigs[:, k],
                out_png,
                title=f"{name}: DESI predicted vs Abacus true",
                labels=("DESI pred", "Abacus true"),
            )
        else:
            # fallback: just DESI histogram
            fig, ax = plt.subplots(1, 1, figsize=(6.5, 4.0))
            ax.hist(preds_lam[:, k], bins=80, density=True, alpha=0.8)
            ax.set_title(f"{name}: DESI predicted")
            ax.set_xlabel("eigenvalue")
            ax.set_ylabel("density")
            fig.tight_layout()
            fig.savefig(out_png, dpi=150)
            plt.close(fig)

    _bar_class_fractions(desi_fracs, out_dir / f"class_fractions_lambda_thr{thr:g}.png", title=f"DESI predicted classes (thr={thr:g})")

    # Also write a small table for downstream joins without FITS.
    np.savez_compressed(
        out_dir / "desi_wedge_index_and_preds.npz",
        global_node_id=desi_global_ids,
        ra=ra,
        dec=dec,
        z=zz,
        lambda1=preds_lam[:, 0],
        lambda2=preds_lam[:, 1],
        lambda3=preds_lam[:, 2],
        cls=desi_cls,
    )

    print(f"Wrote inference outputs to: {out_dir}")


if __name__ == "__main__":
    main()

