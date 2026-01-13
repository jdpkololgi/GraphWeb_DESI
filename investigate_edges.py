import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import torch
import os
from torch_geometric.data import Data

# Paths
# POINTS_PATH = "/pscratch/sd/d/dkololgi/abacus/abacus_cartesian_coords.npy"
# EDGES_PATH = "/pscratch/sd/d/dkololgi/abacus/abacus_delaunay_edges_combined_idx.npy"
CACHE_PATH = "/global/homes/d/dkololgi/GraphWeb_DESI/cache/DESI_alpha_geom.pt"

def main():
    print("Loading data...")
    # try:
    #     points = np.load(POINTS_PATH).astype(np.float64)
    #     edges = np.load(EDGES_PATH)
    # except FileNotFoundError as e:
    #     print(f"Error loading data: {e}")
    #     return

    if os.path.exists(CACHE_PATH):
        print(f"Loading cached graph from {CACHE_PATH}...")
        data = torch.load(CACHE_PATH, weights_only=False)
        points = data.pos.numpy().astype(np.float64)
        edges = data.edge_index.t().numpy()
    else:
        print(f"Cache file not found: {CACHE_PATH}")
        return

    print(f"Loaded {len(points):,} points and {len(edges):,} edges")

    print("Calculating edge lengths...")
    # Calculate lengths in chunks to save memory
    lengths = np.zeros(len(edges), dtype=np.float32)
    chunk_size = 1_000_000
    
    for i in range(0, len(edges), chunk_size):
        end = min(i + chunk_size, len(edges))
        batch_edges = edges[i:end]
        
        # Vectorized distance calculation for the batch
        # points has shape (N, 4), we take first 3 columns for X, Y, Z
        p1 = points[batch_edges[:, 0], :3]
        p2 = points[batch_edges[:, 1], :3]
        
        # Calculate Euclidean distance
        diff = p1 - p2
        dist = np.linalg.norm(diff, axis=1)
        lengths[i:end] = dist
        
        if i % (10 * chunk_size) == 0:
            print(f"Processed {i:,} / {len(edges):,} edges")

    print(f"Processing complete.")
    print(f"Max length: {lengths.max():.2f} Mpc")
    print(f"Min length: {lengths.min():.2f} Mpc")
    print(f"Mean length: {lengths.mean():.2f} Mpc")
    print(f"Median length: {np.median(lengths):.2f} Mpc")

    # Plot 1: Edge Length Distribution
    print("Plotting edge length distribution...")
    plt.figure(figsize=(12, 8))
    plt.hist(lengths, bins=100, log=True, color='skyblue', edgecolor='black', alpha=0.7)
    plt.title("Edge Length Distribution (Log Scale)", fontsize=16)
    plt.xlabel("Length (Mpc)", fontsize=14)
    plt.ylabel("Count", fontsize=14)
    plt.grid(True, alpha=0.3, which="both")
    
    # Add vertical line for 1000 Mpc
    plt.axvline(x=1000, color='red', linestyle='--', label='1000 Mpc Threshold')
    plt.legend()
    
    plt.savefig("edge_length_distribution.png", dpi=300)
    print("Saved edge_length_distribution.png")

    # Plot 2: Spatial Distribution of Long Edges
    long_threshold = 1000 # Mpc
    long_mask = lengths > long_threshold
    long_mask_100 = lengths > 100

    num_edges = len(edges)
    num_edges_100 = len(np.where(long_mask_100)[0])

    long_edges_indices = np.where(long_mask)[0]
    long_edges_indices_100 = np.where(long_mask_100)[0]

    num_long = len(long_edges_indices)
    num_long_100 = len(long_edges_indices_100)  
    
    print(f"Found {num_edges:,} edges")
    print(f"Found {num_edges_100:,} edges longer than 100 Mpc")
    print(f"Found {num_long:,} edges longer than {long_threshold} Mpc ({(num_long/num_edges)*100:.2f}%)")
    print(f"Found {num_long_100:,} edges longer than 100 Mpc ({(num_long_100/len(edges))*100:.2f}%)")

    if num_long > 0:
        # Subsample for visualization
        sample_size = min(10000, num_long)
        print(f"Subsampling {sample_size} long edges for visualization...")
        sample_indices = np.random.choice(long_edges_indices, sample_size, replace=False)
        sample_edges = edges[sample_indices]
        
        print(f"Plotting spatial distribution...")
        
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        
        # Prepare lines for plotting
        lines = []
        for edge in sample_edges:
            p1 = points[edge[0], :3]
            p2 = points[edge[1], :3]
            lines.append([p1, p2])
        
        lc = Line3DCollection(lines, colors='red', linewidths=0.1, alpha=0.3)
        ax.add_collection3d(lc)

        # Set limits based on all points to keep perspective
        ax.set_xlim(points[:, 0].min(), points[:, 0].max())
        ax.set_ylim(points[:, 1].min(), points[:, 1].max())
        ax.set_zlim(points[:, 2].min(), points[:, 2].max())
        
        ax.set_xlabel('X (Mpc)')
        ax.set_ylabel('Y (Mpc)')
        ax.set_zlabel('Z (Mpc)')
        ax.set_title(f"Spatial Distribution of {sample_size} Edges > {long_threshold} Mpc", fontsize=16)
        
        # Add a few points to show the cloud density (very sparse)
        # sparse_points_idx = np.random.choice(len(points), 10000, replace=False)
        # ax.scatter(points[sparse_points_idx, 0], points[sparse_points_idx, 1], points[sparse_points_idx, 2], s=0.1, c='k', alpha=0.1)
        
        plt.savefig("long_edges_spatial.png", dpi=300)
        print("Saved long_edges_spatial.png")
    else:
        print("No long edges found to plot.")

if __name__ == "__main__":
    main()
