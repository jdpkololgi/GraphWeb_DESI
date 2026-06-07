# Deprecated DESI graph / inference stacks

Move artifacts here with `workflows/catalog/deprecate_mpc_h_stack.sh` after the Mpc-parity rebuild validates.

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

Rebuild: `workflows/catalog/rerun_desi_bgs_bright_graph_wedge_mpc.sh`
