# Jraph regression caches (DESI bright wedge + staged mock wedge)

> **Staged mock wedge (15-d SBI cache, graph stack, training):** canonical documentation lives in  
> [`TNG/Illustris/workflows/abacus_tweb/BUILD_JRAPH_CACHES.md`](../../../TNG/Illustris/workflows/abacus_tweb/BUILD_JRAPH_CACHES.md)  
> (Illustris repo). This file covers DESI inference and Abacus wedge calibration paths only.

`TNG/Illustris/workflows/jraph/jraph_pipeline.py` in **`--prediction_mode regression`** loads an SBI pickle via `--cache_path`. Required keys:

| Key | Type | Purpose |
|-----|------|---------|
| `graph` | `jraph.GraphsTuple` | Bidirectional edges (default), 7 node features, 5 edge features |
| `regression_targets` | `jnp.ndarray [N, D]` | Scaled training targets (`D=3` increments or `D=15` enhanced) |
| `eigenvalues_raw` | `np.ndarray [N, 3]` | Physical λ₁, λ₂, λ₃ for metrics / inverse transforms |
| `classification_labels` | `jnp.ndarray [N]` | Discrete CWEB classes 0–3 (void/wall/filament/cluster) |
| `masks` | `(train, val, test)` bool | 70% / 21% / 9% stratified split (default seed 42) |
| `target_scaler` | `sklearn.StandardScaler` | Fitted on **train** split only |
| `stats` | `dict` or `None` | Bounds for 3-d transformed-increment mode |
| `regression_targets_raw` | `np.ndarray [N, D]` | Unscaled assembled targets (Abacus builder only) |

**Discrete classes:** count of λᵢ > 0.2 (same as `build_staged_mock_wedge_truth_npz.classify_eigs`). Stored as `classification_labels`; NPZ truth may use key `cls`.

**15-d vs 3-d:** If the targets FITS has derivative columns (`DLAM*_D*`, `LAP_LAM*`), `build_abacus_sbi_cache.py` builds 15-d targets and sets `use_transformed_eig=false` in training. NPZ-only truth (staged mock) yields **3-d transformed increments** unless you add derivatives.

---

## A. Abacus training wedge (reference — cache exists)

| Artifact | Path |
|----------|------|
| GNN arrays | `/pscratch/sd/d/dkololgi/abacus/graph_constructions/abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_cugraph_gnn_arrays.npz` |
| GNN metadata | `.../abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_cugraph_gnn_metadata.json` |
| **SBI cache** | `/pscratch/sd/d/dkololgi/abacus/sbi_caches/abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_sbi_cache.pkl` |

**Train / eval on Abacus wedge (GPU node):**

```bash
unset PYTHONPATH PYTHONHOME; export PYTHONNOUSERSITE=1
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cosmic_env
cd /global/homes/d/dkololgi/TNG/Illustris

srun --gpus=1 python workflows/jraph/jraph_pipeline.py \
  --prediction_mode regression \
  --cache_path /pscratch/sd/d/dkololgi/abacus/sbi_caches/abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_sbi_cache.pkl \
  --output_dir /pscratch/sd/d/dkololgi/abacus/jraph_runs/wedge_rs7_15d_regression \
  --epochs 4000 --no_transformed_eig
```

(`--no_transformed_eig` matches 15-d Abacus caches.)

**Rebuild cache (CPU, if needed):**

```bash
srun -n1 -c32 python workflows/abacus_tweb/build_abacus_sbi_cache.py \
  --gnn-metadata-path /pscratch/sd/d/dkololgi/abacus/graph_constructions/abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_cugraph_gnn_metadata.json \
  --targets-catalog-path /pscratch/sd/d/dkololgi/abacus/graph_constructions/abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_wedge_targets.fits \
  --output-cache-path /pscratch/sd/d/dkololgi/abacus/sbi_caches/abacus_delaunay_wedge_ra120_140_dec16p5_26p7_z0p25_0p30_rs7_15d_sbi_cache.pkl
```

---

## B. DESI BGS expanded bright wedge (LOA) — inference, not training

**Canonical DESI wedge:** use the expanded Mpc-parity manifest in
`JRAPH_INPUTS_expanded_wedge.txt`. The older narrow bright-wedge manifest
(`JRAPH_INPUTS_bright_wedge.txt`) points at Mpc/h-era artifacts and should be
treated as legacy unless rebuilt with `--coord-units mpc`.

| Artifact | Path |
|----------|------|
| Input manifest | `GraphWeb_DESI/workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt` |
| Full graph dir | `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc` |
| Wedge dir | `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph` |
| GNN arrays | `.../desi_delaunay_wedge_expanded_*_bright_mpc_gnn_arrays.npz` |
| Catalog | `.../desi_delaunay_wedge_expanded_*_bright_mpc_wedge_catalog_minimal.npz` (RA/DEC/Z only — **no truth λ**) |

**No DESI regression cache is required for inference.** LOA has no T-Web eigenvalues; the model predicts them. Use the **Abacus wedge calibration cache** for `target_scaler` and optional truth comparison on the Abacus mock population:

```bash
source /global/homes/d/dkololgi/GraphWeb_DESI/workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt
```

**Graph construction chain:** build the bright catalog with
`workflows/catalog/build_bgs_maglim_catalog.py`, build the full Mpc graph with
`workflows/graph_construction/build_desi_bgs_gudhi_graph.py`, export features
with `workflows/graph_construction/desi_graph_features_cugraph.py`, and induce
the wedge with `workflows/graph_construction/subset_desi_graph_wedge.py`.

On Perlmutter, `workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh` wraps
those three stages as `STEP=2|2b|3` (or `SUBMIT_CHAIN=1` from login). Treat its
hard-coded output dirs / sky cut as the older narrow bright wedge unless you
retarget them to the expanded Mpc paths in `JRAPH_INPUTS_expanded_wedge.txt`
(see `RUNBOOK.md`).

**Inference prerequisites from the code path:**

- The calibration cache should contain the training `target_scaler` and, for the
  validated 15-d run, `node_feature_scaler`.
- DESI edges are duplicated to the bidirectional Abacus cache convention before
  the Jraph forward pass.
- Edge length and density-contrast columns are log-transformed and standardized
  from the Abacus GNN NPZ; pass `--abacus-gnn-arrays`.

**Inference (CPU or GPU):** see
`GraphWeb_DESI/workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`,
`JRAPH_INPUTS_expanded_wedge.txt`, and the thin wrapper
`workflows/catalog/run_infer_expanded_wedge_mpc.sh`.

```bash
conda activate cosmic_env
cd /global/homes/d/dkololgi/GraphWeb_DESI
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

python workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py \
  --abacus-run-dir "${ABACUS_RUN_DIR}" \
  --calibration-cache "${CALIBRATION_CACHE}" \
  --abacus-gnn-arrays "${ABACUS_WEDGE_GNN_ARRAYS}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --output-dir "$(dirname "${INFER_DIR}")" \
  --run-name "${INFER_RUN_NAME}"
```

**Training `jraph_pipeline` on DESI** needs per-node truth λ (cross-match to Abacus annotated CutSky or mock truth). That cache is **not** built here.

**Optional eval cache (future):** cross-match DESI wedge nodes to `staged_mock_wedge_stage3_postcollision_rs7.npz` or annotated CutSky by (RA,DEC,Z), then run `build_abacus_sbi_cache.py` with `--targets-npz-path` on DESI `*_cugraph_gnn_metadata.json`.

---

## C. Staged mock wedge (Illustris — see canonical doc)

Build graph stack, 15-d SBI cache (`*_sbi_cache_15d.pkl`), and training commands in  
[`TNG/Illustris/workflows/abacus_tweb/BUILD_JRAPH_CACHES.md`](../../../TNG/Illustris/workflows/abacus_tweb/BUILD_JRAPH_CACHES.md).

Quick launch from Illustris repo:

```bash
cd /global/homes/d/dkololgi/TNG/Illustris
sbatch workflows/abacus_tweb/submit_staged_mock_wedge_jraph_stack.slurm
```
