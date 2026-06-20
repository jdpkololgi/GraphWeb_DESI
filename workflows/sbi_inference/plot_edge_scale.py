#!/usr/bin/env python3
"""Edge-scale test: is the cluster deficit's fixed-node-feature residual a GRAPH
SCALE shift? DESI has ~12% more galaxies in the wedge -> shorter typical edges. The
GNN's edge scaler (log + StandardScaler) is fit on Abacus, so if DESI edges are
systematically shorter the scaled edge_length is shifted off N(0,1), changing what
the GNN sees even at fixed node features. Compares raw + scaled edge_length.
"""
from __future__ import annotations
import importlib.util, os, sys
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

_ILL = Path(os.environ.get("ILLUSTRIS_ROOT", "/global/homes/d/dkololgi/TNG/Illustris")).resolve()
if str(_ILL) not in sys.path:
    sys.path.insert(0, str(_ILL))
from shared.plot_style import apply_style, ACCENT_COLORS  # noqa: E402

GW = Path("/global/u2/d/dkololgi/GraphWeb_DESI")
spec = importlib.util.spec_from_file_location("agp", GW / "shared" / "abacus_gnn_parity.py")
agp = importlib.util.module_from_spec(spec); spec.loader.exec_module(agp)

PATH1 = "/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz"
DESI = "/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph/desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_gnn_arrays.npz"


def main():
    apply_style()
    a = np.load(PATH1); d = np.load(DESI)
    a_len = a["edge_attr"][:, 0].astype(np.float64); d_len = d["edge_attr"][:, 0].astype(np.float64)
    print(f"galaxies: Abacus {a['x'].shape[0]}  DESI {d['x'].shape[0]}  (DESI/Abacus = {d['x'].shape[0]/a['x'].shape[0]:.3f})")
    print(f"raw edge_length [Mpc]: Abacus median {np.median(a_len):.2f}  DESI median {np.median(d_len):.2f}  (ratio {np.median(d_len)/np.median(a_len):.3f})")

    # scaler fit on Abacus (exactly as inference), applied to both
    sc = agp.fit_edge_length_density_scaler_from_gnn_npz(PATH1, make_bidirectional=True)
    def scaled_len(arr):
        ei = arr["edge_index"]; ea = arr["edge_attr"]
        _, e2 = agp.prepare_edges_for_jraph_forward(ei, ea, sc, make_bidirectional=True)
        return e2[:, 0]
    a_s = scaled_len(a); d_s = scaled_len(d)
    print(f"scaled log(edge_length): Abacus mean {a_s.mean():.3f} std {a_s.std():.3f}  |  DESI mean {d_s.mean():.3f} std {d_s.std():.3f}")
    print(f"  => DESI edge-scale offset from training: {d_s.mean():.3f} sigma (negative = shorter/denser than Abacus)")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    bins = np.linspace(0, np.percentile(np.concatenate([a_len, d_len]), 99), 80)
    axes[0].hist(a_len, bins=bins, density=True, histtype="step", lw=2, color="#9a9a93", label=f"Abacus (med {np.median(a_len):.1f})")
    axes[0].hist(d_len, bins=bins, density=True, histtype="step", lw=2, color=ACCENT_COLORS["magenta"], label=f"DESI (med {np.median(d_len):.1f})")
    axes[0].set_xlabel("raw edge length [Mpc]"); axes[0].set_ylabel("density"); axes[0].set_title("Graph edge lengths"); axes[0].legend()
    bb = np.linspace(-4, 4, 80)
    axes[1].hist(a_s, bins=bb, density=True, histtype="step", lw=2, color="#9a9a93", label=f"Abacus (mean {a_s.mean():.2f})")
    axes[1].hist(d_s, bins=bb, density=True, histtype="step", lw=2, color=ACCENT_COLORS["magenta"], label=f"DESI (mean {d_s.mean():.2f})")
    axes[1].axvline(0, color="#888", ls=":", lw=1)
    axes[1].set_xlabel("scaled log(edge length)  [Abacus-fit scaler]"); axes[1].set_ylabel("density")
    axes[1].set_title("What the GNN sees (Abacus=N(0,1) by construction)"); axes[1].legend()
    fig.suptitle("Edge-scale domain shift: DESI graph is denser -> edges scaled off training")
    out = "/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear/edge_scale.png"
    fig.savefig(out, bbox_inches="tight", dpi=200); plt.close(fig); print("Saved:", out)


if __name__ == "__main__":
    main()
