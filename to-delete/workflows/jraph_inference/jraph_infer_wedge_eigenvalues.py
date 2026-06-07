#!/usr/bin/env python3
"""Run Abacus-trained Jraph regression model on a DESI wedge inference cache.

Inputs:
- --input-cache: output of workflows/jraph_inference/build_desi_wedge_jraph_cache.py
- --model-path: trained Jraph model params/checkpoint (.pkl)
- --calibration-cache: Abacus wedge SBI cache used for training (provides target_scaler)

Outputs:
- predictions pickle with galaxy_table + predicted eigenvalues (λ1,λ2,λ3)
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

# Default to CPU-safe; override with CUDA_VISIBLE_DEVICES on GPU nodes.
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")

# Allow canonical workflow scripts to resolve repo-root modules after reorganization.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import jax
import haiku as hk

from shared.graph_net_models import make_graph_network


def _load_params(model_path: Path) -> tuple[object, int | None]:
    with model_path.open("rb") as f:
        obj = pickle.load(f)
    if isinstance(obj, dict) and "params" in obj:
        return obj["params"], obj.get("epoch")
    return obj, None


def _eig_from_15(raw15: np.ndarray) -> np.ndarray:
    """Reconstruct (λ1, λ2, λ3) from first three raw channels (λ1, Δλ2, Δλ3)."""
    raw15 = np.asarray(raw15, dtype=np.float64)
    l1 = raw15[:, 0]
    d2 = np.maximum(raw15[:, 1], 1e-7)
    d3 = np.maximum(raw15[:, 2], 1e-7)
    l2 = l1 + d2
    l3 = l2 + d3
    return np.stack([l1, l2, l3], axis=-1)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-cache", type=Path, required=True)
    p.add_argument("--model-path", type=Path, required=True)
    p.add_argument("--calibration-cache", type=Path, required=True)
    p.add_argument("--output-path", type=Path, required=True)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--latent_size", type=int, default=96)
    p.add_argument("--num_heads", type=int, default=8)
    p.add_argument("--num_passes", type=int, default=8)
    p.add_argument("--dropout", type=float, default=0.15)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")

    with args.input_cache.expanduser().resolve().open("rb") as f:
        in_cache = pickle.load(f)
    graph = in_cache["graph"]
    gal = in_cache["galaxy_table"]

    with args.calibration_cache.expanduser().resolve().open("rb") as f:
        cal = pickle.load(f)
    target_scaler = cal.get("target_scaler")
    if target_scaler is None:
        raise KeyError("Calibration cache missing target_scaler; cannot inverse-transform predictions.")

    params, ckpt_epoch = _load_params(args.model_path.expanduser().resolve())

    net_fn = make_graph_network(
        num_passes=args.num_passes,
        latent_size=args.latent_size,
        num_heads=args.num_heads,
        dropout_rate=args.dropout,
        output_dim=15,
    )
    net = hk.transform(net_fn)

    @jax.jit
    def predict(p, g, rng):
        return net.apply(p, rng, g, is_training=False).nodes

    rng = jax.random.PRNGKey(args.seed + 999)
    preds_scaled = np.asarray(predict(params, graph, rng))
    if preds_scaled.ndim != 2 or preds_scaled.shape[1] != 15:
        raise ValueError(f"Unexpected preds shape {preds_scaled.shape}; expected [N,15].")

    preds_raw15 = target_scaler.inverse_transform(preds_scaled)
    preds_lambda = _eig_from_15(preds_raw15)

    out = args.output_path.expanduser().resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "timestamp": ts,
        "model_path": str(args.model_path),
        "checkpoint_epoch": ckpt_epoch,
        "input_cache": str(args.input_cache),
        "calibration_cache": str(args.calibration_cache),
        "hparams": {
            "latent_size": args.latent_size,
            "num_heads": args.num_heads,
            "num_passes": args.num_passes,
            "dropout": args.dropout,
            "seed": args.seed,
        },
        "preds_scaled_15d": preds_scaled,
        "preds_raw_15d": np.asarray(preds_raw15, dtype=np.float64),
        "preds_lambda123": np.asarray(preds_lambda, dtype=np.float64),
        "galaxy_table": gal,
        "wedge_bounds": in_cache.get("wedge_bounds"),
    }
    with out.open("wb") as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)

    meta = {
        "timestamp": ts,
        "model_path": str(args.model_path),
        "checkpoint_epoch": ckpt_epoch,
        "input_cache": str(args.input_cache),
        "calibration_cache": str(args.calibration_cache),
        "output_path": str(out),
        "n_galaxies": int(getattr(gal, "shape", [0])[0]),
    }
    meta_path = out.with_suffix(".meta.json")
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print(f"Wrote predictions: {out}")
    print(f"Wrote meta: {meta_path}")


if __name__ == "__main__":
    main()

