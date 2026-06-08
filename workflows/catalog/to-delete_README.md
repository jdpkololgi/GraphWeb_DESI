# Deprecated DESI graph / inference stacks

This folder is for retired DESI wedge artifacts after the Mpc-parity rebuild
validates.

## Why deprecated

| Issue | Effect |
|-------|--------|
| Coordinates in **Mpc/h** (`× h`) vs Abacus **Mpc** | Mean edge lengths ~0.71×; density / inertia skewed |
| Inference without **log + StandardScaler** on edge cols 0,4 | ~71% void CWEB, ~71% negative Δλ₃ |
| Pre-`eigfix` increment mislabel | ~87% false walls |

## Replacement paths

- Full graph: `outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc`
- Expanded wedge: `outputs/desi_wedge_expanded_*_bright_mpc_from_fullgraph`
- Inference: `inference_outputs/infer_expanded_*_eigfix_edgescale_mpc`

Rebuild with the Python workflow chain:

1. `workflows/catalog/build_bgs_maglim_catalog.py`
2. `workflows/graph_construction/build_desi_bgs_gudhi_graph.py --coord-units mpc`
3. `workflows/graph_construction/desi_graph_features_cugraph.py`
4. `workflows/graph_construction/subset_desi_graph_wedge.py`
5. `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`
