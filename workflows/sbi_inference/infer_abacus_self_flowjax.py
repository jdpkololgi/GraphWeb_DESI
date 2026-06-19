#!/usr/bin/env python3
"""Run the linear FlowJAX NPE on its OWN Abacus wedge (self-inference) to get the
'Abacus NPE' reference: per-galaxy posterior-mean λ, embeddings, and a posterior
sample subset on the Abacus TEST split — parallel arrays to the DESI inference, so
the 3-way eigenvalue distributions (truth vs Abacus-NPE vs DESI-NPE) and the
embedding-space comparison can be built. The Abacus graph is already preprocessed
in the cache, so no node/edge transforms are needed here.
"""
from __future__ import annotations
import argparse, os, pickle, sys
from pathlib import Path
import numpy as np

os.environ.setdefault("PYTHONNOUSERSITE", "1")
os.environ.setdefault("XLA_PYTHON_CLIENT_PREALLOCATE", "false")
os.environ.setdefault("XLA_PYTHON_CLIENT_ALLOCATOR", "platform")
_ILLUSTRIS_ROOT = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILLUSTRIS_ROOT) not in sys.path:
    sys.path.insert(0, str(_ILLUSTRIS_ROOT))
import jax
from workflows.sbi.plot_flowjax_posteriors import load_flowjax_model, create_gnn_and_flow, batched_sample_posterior
from shared.eigenvalue_transformations import samples_to_raw_eigenvalues, posterior_to_classprobs

CLASS_ORDER = ("void", "wall", "filament", "cluster")


def main(args):
    out_dir = args.output_dir.expanduser().resolve(); out_dir.mkdir(parents=True, exist_ok=True)
    print("JAX devices:", jax.devices(), flush=True)
    gnn_params, config, target_scaler, flow_filename, increment_mode = load_flowjax_model(str(args.model_path))
    with args.calibration_cache.open("rb") as f:
        cache = pickle.load(f)
    graph = cache["graph"]
    eig_raw = np.asarray(cache["eigenvalues_raw"], dtype=np.float64)
    test_mask = np.asarray(cache["masks"][2])
    idx = np.where(test_mask)[0]

    gnn, flow = create_gnn_and_flow(config, flow_filename, graph, jax.random.key(args.seed))
    emb_all = np.asarray(gnn.apply(gnn_params, jax.random.key(0), graph, is_training=False))
    emb = emb_all[idx]
    print(f"Abacus test nodes: {len(idx)}; embeddings {emb.shape}", flush=True)

    samples = batched_sample_posterior(flow, emb, args.num_posterior_samples, jax.random.key(123), chunk_size=512)
    samples_raw = samples_to_raw_eigenvalues(samples, target_scaler, increment_mode)
    samples_raw = np.sort(samples_raw, axis=-1)  # post-hoc ascending (matches DESI run)
    lambda_mean = samples_raw.mean(axis=1)
    cp = posterior_to_classprobs(samples_raw, lambda_th=args.lambda_threshold)
    classprob = np.stack([np.asarray(cp[c]) for c in CLASS_ORDER], axis=1)
    hard = classprob.argmax(axis=1).astype(np.int8)

    rng = np.random.default_rng(0)
    sub = np.sort(rng.choice(len(idx), min(len(idx), args.save_sample_subset), replace=False))
    np.savez_compressed(
        out_dir / "abacus_self_flowjax_preds.npz",
        node_index=idx, eig_truth=eig_raw[idx], lambda_mean=lambda_mean,
        classprob=classprob, hard_class=hard, embeddings=emb.astype(np.float32),
        sample_subset_idx=sub, lambda_samples_subset=samples_raw[sub].astype(np.float32))
    print("DESI-comparison ready. Abacus NPE class fracs:",
          {CLASS_ORDER[k]: round(float(np.mean(hard == k)), 3) for k in range(4)}, flush=True)
    print("Saved:", out_dir / "abacus_self_flowjax_preds.npz", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-path", type=Path, required=True)
    ap.add_argument("--calibration-cache", type=Path, required=True)
    ap.add_argument("--num-posterior-samples", type=int, default=128)
    ap.add_argument("--lambda-threshold", type=float, default=0.2)
    ap.add_argument("--save-sample-subset", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--output-dir", type=Path,
                    default=Path("/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear"))
    main(ap.parse_args())
