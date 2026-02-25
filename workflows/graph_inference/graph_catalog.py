import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

# Fast-exit help path so `--help` works even if heavy deps are unavailable.
if __name__ == "__main__" and any(arg in ("-h", "--help") for arg in sys.argv[1:]):
    print("Usage: python graph_catalog.py [--graph-type {alpha,delaunay}] [--cache-mode {rebuild,prefer-cache,cache-only}] [--no-summary-plot] [--cache-dir PATH] [--model-path PATH] [--scaler-path PATH] [--vac-output-path PATH] [--reference-catalog-path PATH]")
    raise SystemExit(0)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scienceplots
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.preprocessing import PowerTransformer
from torch_geometric.nn import GATv2Conv
from torch_geometric.utils import from_networkx

# Allow canonical workflow scripts to resolve repo-root modules after reorganization.
REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from config_paths import (
    GRAPHWEB_CACHE_DIR,
    GRAPHWEB_VAC_OUTPUT_PATH,
    ILLUSTRIS_GAT_MODEL_PATH,
    ILLUSTRIS_REPO_ROOT,
    ILLUSTRIS_SCALER_PATH,
    TNG_REFERENCE_CATALOG_PATH,
)

# Workflow status: ACTIVE (canonical GraphWeb DESI inference pipeline)

sys.path.append("../")
if ILLUSTRIS_REPO_ROOT not in sys.path:
    sys.path.append(ILLUSTRIS_REPO_ROOT)
from Network_stats import network
from Utilities import cat


plt.style.use(["science", "no-latex", "dark_background"])


@dataclass(frozen=True)
class CachePaths:
    graph: str
    geom: str
    features: str
    zcat: str


def build_cache_paths(cache_dir: str, graph_type: str) -> CachePaths:
    if graph_type == "alpha":
        return CachePaths(
            graph=os.path.join(cache_dir, "DESI_alpha_graph.pt"),
            geom=os.path.join(cache_dir, "DESI_alpha_geom.pt"),
            features=os.path.join(cache_dir, "DESI_alpha_features.pt"),
            zcat=os.path.join(cache_dir, "DESI_NETWORKalpha_zcat.pt"),
        )
    return CachePaths(
        graph=os.path.join(cache_dir, "DESI_delaunay_graph.pt"),
        geom=os.path.join(cache_dir, "DESI_delaunay_geom.pt"),
        features=os.path.join(cache_dir, "DESI_delaunay_features.pt"),
        zcat=os.path.join(cache_dir, "DESI_NETWORK_delaunay_zcat.pt"),
    )


def cache_exists(paths: CachePaths) -> bool:
    return (
        os.path.exists(paths.graph)
        and os.path.exists(paths.geom)
        and os.path.exists(paths.features)
        and os.path.exists(paths.zcat)
    )


def load_cached_bundle(paths: CachePaths):
    print("Loading cached data...")
    graph = torch.load(paths.graph, weights_only=False)
    desi_geom = torch.load(paths.geom, weights_only=False)
    desi_features = pd.read_pickle(paths.features)
    zcat = pd.read_pickle(paths.zcat)
    print("Cached data loaded successfully.")
    return graph, desi_geom, desi_features, zcat


def save_with_memory_check(obj, path, obj_name):
    import gc
    import pickle
    import psutil

    try:
        available_memory = psutil.virtual_memory().available / (1024**3)
        print(f"Available memory before saving {obj_name}: {available_memory:.2f} GB")
        if available_memory < 2.0:
            print(f"Warning: Low memory before saving {obj_name}")
            gc.collect()

        print(f"Saving {obj_name} to {path}...")
        torch.save(obj, path)
        print(f"Successfully saved {obj_name}")
        del obj
        gc.collect()
    except Exception as e:
        print(f"Error saving {obj_name}: {e}")
        try:
            with open(path, "wb") as f:
                pickle.dump(obj, f, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Successfully saved {obj_name} using pickle")
        except Exception as e2:
            print(f"Failed to save {obj_name} with both methods: {e2}")
            return False
    return True


def build_desi_graph_bundle(graph_type: str, first_moment_matching: bool, scaler_path: str):
    print("Cached data missing or incomplete, creating new objects...")
    print(f"Creating {graph_type} network object for DESI BGS galaxies...")
    desi_network = network(masscut=9.0, from_DESI=True)
    print("DESI network object created")
    zcat = desi_network.DESI_GAL_CAT.zcat.to_pandas()

    if graph_type == "alpha":
        graph = desi_network.galaxy_alpha_complex_network(xyzplot=False)
        print("DESI alpha complex graph created")
        desi_network.network_stats_alpha(G=graph)
        print("DESI alpha complex network stats calculated")
    else:
        graph = desi_network.subhalo_delaunay_network(xyzplot=False)
        print("DESI delaunay graph created")
        desi_network.network_stats_delaunay()
        print("DESI delaunay network stats calculated")

    desi_geom = from_networkx(graph, group_edge_attrs=["length"])
    print("DESI graph converted to torch_geometric Data object")

    if first_moment_matching:
        scaler = torch.load(scaler_path, weights_only=False)
        features_data = scaler.transform(desi_network.data + 1e-6)
    else:
        scaler = PowerTransformer(method="box-cox")
        features_data = scaler.fit_transform(desi_network.data + 1e-6)

    desi_features = pd.DataFrame(
        features_data,
        index=desi_network.data.index,
        columns=desi_network.data.columns,
    )
    print("Applying Domain Adaptation: Centering DESI features to Mean=0...")
    desi_features = desi_features - desi_features.mean()
    desi_geom.x = torch.tensor(desi_features.values, dtype=torch.float32)
    print("DESI features scaled and converted to torch tensor")
    return graph, desi_geom, desi_features, zcat, desi_network


class SimpleGAT(nn.Module):
    def __init__(self, input_dim, output_dim, num_heads=1):
        super().__init__()
        hidden_dim = 20
        total_hidden = hidden_dim * num_heads
        self.gat_layer1 = GATv2Conv(input_dim, hidden_dim, edge_dim=1, heads=num_heads, concat=True)
        self.gat_layer2 = GATv2Conv(total_hidden, hidden_dim, edge_dim=1, heads=num_heads, concat=True)
        self.gat_layer3 = GATv2Conv(total_hidden, hidden_dim, edge_dim=1, heads=num_heads, concat=True)
        self.gat_layer4 = GATv2Conv(hidden_dim * num_heads, output_dim, edge_dim=1, heads=1, concat=False)
        self.gat_dropout = nn.Dropout(p=0.1)
        self.norm1 = nn.LayerNorm(total_hidden)
        self.norm2 = nn.LayerNorm(total_hidden)
        self.norm3 = nn.LayerNorm(total_hidden)

    def forward(self, x, edge_index, edge_weight=None):
        x1 = self.gat_dropout(F.relu(self.norm1(self.gat_layer1(x, edge_index, edge_attr=edge_weight))))
        x2 = self.gat_dropout(F.relu(self.norm2(self.gat_layer2(x1, edge_index, edge_attr=edge_weight)) + x1))
        x3 = self.gat_dropout(F.relu(self.norm3(self.gat_layer3(x2, edge_index, edge_attr=edge_weight)) + x2))
        return self.gat_layer4(x3, edge_index, edge_attr=edge_weight)


def run_inference(desi_geom, model_path: str):
    model = SimpleGAT(10, 4, num_heads=4)
    model.load_state_dict(torch.load(model_path, map_location="cpu"))
    model.eval()
    print("Inference model loaded and set to eval mode")
    with torch.no_grad():
        out = model(desi_geom.x, desi_geom.edge_index, edge_weight=desi_geom.edge_attr)
        preds = out.argmax(dim=1).numpy()
        probs = F.softmax(out, dim=1).numpy()
    print("Inference completed, predictions and probabilities obtained")
    return preds, probs


def write_vac(zcat, preds, probs, vac_output_path: str):
    print("Saving pre-release VAC of DESI BGS galaxies with cosmic web predictions...")
    zcat["GAT_ENV"] = preds
    zcat["GAT_VOID_PROB"] = probs[:, 0]
    zcat["GAT_WALL_PROB"] = probs[:, 1]
    zcat["GAT_FILAMENT_PROB"] = probs[:, 2]
    zcat["GAT_CLUSTER_PROB"] = probs[:, 3]
    if hasattr(zcat, "to_pandas"):
        zcat = zcat.to_pandas()
    os.makedirs(os.path.dirname(vac_output_path), exist_ok=True)
    zcat.to_pickle(vac_output_path)
    print(f"Wrote VAC: {vac_output_path}")
    return zcat


def maybe_plot_summary(testcat, preds):
    environ_dicts = {0: "Void", 1: "Wall", 2: "Filament", 3: "Cluster"}
    custom_palette = {0: "#80ffdb", 1: "#3a86ff", 2: "#ff006e", 3: "#ffbe0b"}
    testcat.cweb_classify(xyzplot=False)
    bins = np.arange(5) - 0.5
    desi_hist, _ = np.histogram(preds, bins=bins, density=True)
    tweb_hist, _ = np.histogram(testcat.cweb, bins=bins, density=True)
    bar_width = 0.4
    x = np.arange(4)
    plt.figure(figsize=(10, 6))
    plt.bar(x - bar_width / 2, desi_hist, width=bar_width, color="#80ffdb", edgecolor="black", label="BGS")
    plt.bar(x + bar_width / 2, tweb_hist, width=bar_width, color="#3a86ff", edgecolor="black", alpha=0.7, label="IllustrisTNG T-WEB")
    for i in range(4):
        plt.text(x[i] - bar_width / 2, desi_hist[i] + 0.005, f"{desi_hist[i]*100:.1f}%", ha="center", va="bottom", fontsize=10)
        plt.text(x[i] + bar_width / 2, tweb_hist[i] + 0.005, f"{tweb_hist[i]*100:.1f}%", ha="center", va="bottom", fontsize=10)
    plt.xticks(x, [environ_dicts[i] for i in range(4)])
    plt.xlabel("Cosmic Web Environment")
    plt.ylabel("Frequency")
    plt.title("Distribution of Cosmic Web Environments in DESI Network")
    plt.legend()
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.show()


def parse_args():
    parser = argparse.ArgumentParser(description="Build/load DESI graph and run GAT environment inference.")
    parser.add_argument("--graph-type", choices=["alpha", "delaunay"], default="alpha")
    parser.add_argument(
        "--cache-mode",
        choices=["rebuild", "prefer-cache", "cache-only"],
        default="rebuild",
        help="rebuild: force new graph/cache; prefer-cache: load cache if available; cache-only: fail if cache missing.",
    )
    parser.add_argument("--first-moment-matching", action="store_true")
    parser.add_argument("--cache-dir", default=GRAPHWEB_CACHE_DIR)
    parser.add_argument("--model-path", default=ILLUSTRIS_GAT_MODEL_PATH)
    parser.add_argument("--scaler-path", default=ILLUSTRIS_SCALER_PATH)
    parser.add_argument("--vac-output-path", default=GRAPHWEB_VAC_OUTPUT_PATH)
    parser.add_argument("--no-summary-plot", action="store_true")
    parser.add_argument("--reference-catalog-path", default=TNG_REFERENCE_CATALOG_PATH)
    return parser.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.cache_dir, exist_ok=True)
    cache_paths = build_cache_paths(args.cache_dir, args.graph_type)
    testcat = cat(path=args.reference_catalog_path, snapno=99, masscut=1e9)
    print("Created cat object for TNG300 galaxies")

    can_use_cache = cache_exists(cache_paths)
    if args.cache_mode == "cache-only" and not can_use_cache:
        raise FileNotFoundError("Cache mode is cache-only but required cache files are missing.")

    if args.cache_mode in ("prefer-cache", "cache-only") and can_use_cache:
        graph, desi_geom, desi_features, zcat = load_cached_bundle(cache_paths)
        desi_network = None
    else:
        graph, desi_geom, desi_features, zcat, desi_network = build_desi_graph_bundle(
            graph_type=args.graph_type,
            first_moment_matching=args.first_moment_matching,
            scaler_path=args.scaler_path,
        )

        print("Saving data to cache...")
        save_with_memory_check(graph, cache_paths.graph, "Graph")
        save_with_memory_check(desi_geom, cache_paths.geom, "DESI_geom")
        try:
            desi_features.to_pickle(cache_paths.features)
            print("Successfully saved DESI_features")
        except Exception as e:
            print(f"Error saving DESI_features: {e}")
        if desi_network is not None:
            try:
                desi_network.DESI_GAL_CAT.zcat.to_pandas().to_pickle(cache_paths.zcat)
                print("Successfully saved DESI_NETWORK.zcat")
            except Exception as e:
                print(f"Error saving DESI_NETWORK.zcat: {e}")
        print("Data cached successfully.")

    preds, probs = run_inference(desi_geom, model_path=args.model_path)
    zcat = write_vac(zcat, preds, probs, vac_output_path=args.vac_output_path)
    if not args.no_summary_plot:
        maybe_plot_summary(testcat, preds)


if __name__ == "__main__":
    main()
