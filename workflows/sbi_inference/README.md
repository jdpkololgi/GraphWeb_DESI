# FlowJAX/SBI DESI wedge inference

This directory contains the DESI-side posterior inference and diagnostics for the
Abacus-trained FlowJAX neural posterior estimator (NPE). It is the current
DESI wedge path for per-galaxy ordered tidal-tensor eigenvalue posteriors:

```text
DESI Mpc-parity wedge graph
  -> Abacus-parity node/edge transforms
  -> Illustris GNN encoder + FlowJAX posterior
  -> per-galaxy lambda samples, class probabilities, and closure figures
```

This is separate from:

- `workflows/graph_inference/graph_catalog.py`: all-sky GAT VAC-style class
  prediction.
- `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`: Jraph
  point-regression for eigenvalues on a DESI wedge.

## Scientific intent

The model infers the ordered T-web eigenvalues
`lambda_1 <= lambda_2 <= lambda_3`. Void/wall/filament/cluster classes are
derived from posterior samples using `lambda_threshold` (default `0.2`), not
treated as the fundamental target. DESI has no per-galaxy truth labels, so the
core DESI-side checks are truth-free:

- class-probability count/marginal consistency;
- posterior width and entropy maps;
- Abacus-vs-DESI eigenvalue and embedding comparisons;
- galaxy-property closure tests after joining FastSpecFit properties.

Calibration tests that require truth (TARP/SBC and Abacus test-set metrics)
remain in the sibling Illustris training workflow.

## Required inputs and parity constraints

`infer_desi_wedge_flowjax.py` requires a DESI wedge built from the Gudhi/cuGraph
Mpc-parity graph products plus the matching Abacus training artifacts.

Required inputs:

- `--model-path`: FlowJAX model pickle from the Illustris SBI run.
- `--calibration-cache`: training cache with `node_feature_scaler`,
  `eigenvalues_raw`, and target scaler metadata.
- `--abacus-gnn-arrays`: raw path1 fiberassign Abacus wedge GNN arrays used to
  fit the DESI edge scaler.
- `--desi-gnn-arrays`: DESI wedge GNN arrays.
- `--desi-gnn-metadata`: DESI wedge metadata JSON.
- `--desi-global-node-ids`: row indices into the source BGS graph catalog.
- `--desi-wedge-catalog-npz`: minimal RA/Dec/z wedge catalog.

Important invariants enforced by the script:

- DESI coordinate units must be comoving `mpc`; `mpc/h` or unknown metadata
  aborts inference.
- Node features must be exactly:
  `Degree, Clustering, Density, Neigh Density, I_eig1, I_eig2, I_eig3`.
- Edge features are duplicated bidirectionally, then edge length and density
  contrast are log-transformed and standardized with an Abacus path1 wedge
  scaler.
- `ILLUSTRIS_ROOT` must point to the sibling Illustris repo so imports such as
  `shared.graph_net_models` and `workflows.sbi.plot_flowjax_posteriors` resolve
  to the training code.

Use the expanded Mpc wedge manifest as the source of canonical DESI paths:

```bash
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt
```

## Inference

The checked-in launcher captures the validated **baseline** Perlmutter pattern
(no `--scale-invariant-features`):

```bash
bash workflows/sbi_inference/run_infer_desi_wedge_flowjax.sh
```

That launcher still writes `--run-name desi_wedge_flowjax_linear`. Downstream
examples in this README often point at the scale-invariant production directory
`desi_wedge_flowjax_linear_si`; those require a matching SI model/cache and an
explicit `--scale-invariant-features` flag — do not mix the two.

The underlying command shape is:

```bash
export ILLUSTRIS_ROOT=/global/homes/d/dkololgi/TNG/Illustris
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

FLOWJAX_RUN=/pscratch/sd/d/dkololgi/abacus/sbi_runs/path1_wedge_flowjax_3d_Bcorrected_linear
MODEL="${FLOWJAX_RUN}/flowjax_sbi_model_seed_42_YYYYMMDD_HHMMSS.pkl"
CACHE=/pscratch/sd/d/dkololgi/abacus/sbi_caches/path1_flowjax_3d_lineareig/processed_jraph_data_mc1e+09_v2_scaled_3_linear_eig.pkl
PATH1_EDGE_ARRAYS=/pscratch/sd/d/dkololgi/abacus/graph_constructions/wedges/path1_fiberassign/path1_fiberassign_mock_bgs_maglim_rs7_wedge_ra120_160_dec14p5_30p6_z0p2_0p3_cugraph_gnn_arrays.npz

python -u workflows/sbi_inference/infer_desi_wedge_flowjax.py \
  --model-path "${MODEL}" \
  --calibration-cache "${CACHE}" \
  --abacus-gnn-arrays "${PATH1_EDGE_ARRAYS}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --num-posterior-samples 128 \
  --lambda-threshold 0.2 \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs \
  --run-name desi_wedge_flowjax_linear
```

For the scale-invariant production run, replace `MODEL` and `CACHE` with the
matching scale-invariant training artifacts, set
`--run-name desi_wedge_flowjax_linear_si`, and add
`--scale-invariant-features`.

Outputs are written to `--output-dir/--run-name/`:

- `desi_wedge_flowjax_preds.npz`: `global_node_id`, sky coordinates,
  `lambda_mean`, `lambda_std`, `classprob`, `hard_class`, `p_exceed`,
  embeddings, and a posterior-sample subset.
- `summary.json`: provenance, scaler statistics, class fractions, ordering
  checks, and Abacus reference fractions.

### Domain-adaptation and production flags

The inference script exposes three transfer-diagnostic modes:

- `--edge-domain-adapt`: re-standardizes DESI scaled edge length and density
  contrast to training-like `N(0,1)`.
- `--node-domain-adapt`: re-standardizes DESI box-cox node features to
  training-like `N(0,1)`.
- `--scale-invariant-features`: uses per-graph-median normalized node and edge
  scale features. Use this only with a matching scale-invariant training cache
  and model.

Treat node/edge domain adaptation as diagnostics unless the science log records
the corresponding run as the selected data product. The scale-invariant mode is
the durable Route A transfer path, but only when model, cache, and DESI transform
all match.

## Abacus self-inference reference

For DESI comparison plots, run the trained NPE on its own Abacus test wedge:

```bash
python workflows/sbi_inference/infer_abacus_self_flowjax.py \
  --model-path "${MODEL}" \
  --calibration-cache "${CACHE}" \
  --output-dir /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear_si
```

This writes `abacus_self_flowjax_preds.npz`, which supplies Abacus NPE
posterior means, embeddings, and class probabilities for three-way
Abacus-truth / Abacus-NPE / DESI-NPE figures.

## Figures and truth-free diagnostics

After DESI inference, generate the main diagnostic figures on a CPU node:

```bash
RUN_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear_si
ABACUS_SELF=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear_si/abacus_self_flowjax_preds.npz

python workflows/sbi_inference/plot_desi_wedge_flowjax.py \
  --preds-npz "${RUN_DIR}/desi_wedge_flowjax_preds.npz" \
  --summary-json "${RUN_DIR}/summary.json" \
  --calibration-cache "${CACHE}" \
  --abacus-self-npz "${ABACUS_SELF}" \
  --output-dir "${RUN_DIR}" \
  --no-regression
```

Useful plot flags:

- `--no-regression`: drop the off-message "Regression DESI" bars from
  `class_fractions_comparison.png` (SBI-room three-series layout).
- `--only-class-fractions`: regenerate only that bar chart and exit.

Common products include:

- `class_fractions_comparison.png`;
- `class_sky_map.png`;
- `posterior_width_sky_map.png`;
- `eigenvalue_distributions_3way.png` (the on-disk three-way λ overlay; some
  comments/print labels still say `domain_shift_overlay.png`);
- `embedding_pca.png` and optional `embedding_pca_3d.html`;
- optional UMAP and 3D HTML figures when dependencies are available.

For the cluster-deficit transfer story, see:

- `CLUSTER_DEFICIT_FIX_PLAN.md`;
- `plot_edge_scale.py`;
- `plot_cluster_recovery_bars.py`;
- `plot_eigenvalue_corner.py`;
- `plot_eig_dist_buildup.py`.

## Luminosity support diagnostic (experimental)

`desi_absmag_kcorr.py` is a separate DESI-environment diagnostic that computes
official LSS k+e-corrected `ABSMAG_RP1` and prints a marginal comparison with
Abacus `R_MAG_ABS`. It is **not** a production train/inference gate. Launch
commands, runtime constraints, and interpretation limits are in `RUNBOOK.md`.

## Galaxy-property closure join

Join DESI wedge posteriors back to FastSpecFit properties:

```bash
python workflows/sbi_inference/build_desi_wedge_property_join.py \
  --preds "${RUN_DIR}/desi_wedge_flowjax_preds.npz"
```

The join recovers `TARGETID` through `global_node_id` as an index into the source
BGS catalog, averages duplicate hemisphere copies, and writes
`desi_wedge_env_props.parquet` plus `join_report.json`. The report includes
FastSpecFit match fraction, duplicate-target count, valid-sSFR fraction, and an
RA/Dec cross-check.

Then make closure figures:

```bash
python workflows/sbi_inference/plot_property_environment_closure.py \
  --table "${RUN_DIR}/desi_wedge_env_props.parquet" \
  --outdir "${RUN_DIR}"
```

Expected qualitative closure: quenched fraction and rest-frame `g-r` increase
toward denser inferred environments, while median log sSFR decreases.

### CIGALE property swap (Approach A)

If FastSpecFit SFRs look pathological for science plots, re-join CIGALE-HZ
stellar mass and SFR onto the **existing** SI wedge posteriors without
re-running inference:

```bash
python workflows/sbi_inference/build_cigale_rejoin.py
```

Hard-coded defaults (edit the script if paths move):

- input: SI FastSpecFit join parquet under `desi_wedge_flowjax_linear_si/`
- CIGALE catalogue on CFS (`CG_15` fiducial; `CG_5` kept for provenance)
- output: `.../desi_wedge_cigale_hz/desi_wedge_env_props.parquet`

The output keeps the same column schema as the FastSpecFit join so downstream
plot scripts can be repointed via `WEDGE_PARQUET` / `WEDGE_FIGDIR`. Coverage is
~90% of wedge galaxies in the goodPhoto CIGALE catalogue. Environment
inference (graph, posteriors, DELTACHI2>=25 parity) is unchanged.

### Property-science figures

These scripts are CPU-only, take no CLI flags, and read/write via env vars
(defaults shown):

```bash
export WEDGE_PARQUET=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_cigale_hz/desi_wedge_env_props.parquet
export WEDGE_FIGDIR=/pscratch/sd/d/dkololgi/graphweb_desi/figures/desi_wedge_cigale_hz
export ILLUSTRIS_ROOT=/global/homes/d/dkololgi/TNG/Illustris

python workflows/sbi_inference/plot_sfms_environment.py
python workflows/sbi_inference/plot_mstar_color_environment.py
python workflows/sbi_inference/plot_env_mass_continuous.py
python workflows/sbi_inference/investigate_env_property_signal.py
```

Intent:

- `plot_sfms_environment.py`: SFR–M* hexbins by inferred class and continuous
  λ colouring (MS / green-valley / red-sequence bands).
- `plot_mstar_color_environment.py`: M* vs `(g-r)` KDE contours by class.
- `plot_env_mass_continuous.py`: mass-controlled quenched / sSFR / colour
  surfaces and an animated mass slider (defaults to the CIGALE-HZ parquet).
- `investigate_env_property_signal.py`: decompose weak env↔property signal into
  inference noise vs local-density physics (kNN density contrast, mass tertiles,
  Abacus self-eval accuracy, posterior-width checks). Hard-coded to the
  CIGALE-HZ parquet + SI Abacus self-eval npz.

Note: `plot_env_mass_continuous.py` currently inserts a hard-coded Illustris
path for `shared.plot_style` rather than `ILLUSTRIS_ROOT`; set that path if you
run off the usual NERSC home layout. Prefer `ILLUSTRIS_ROOT` for the other
scripts.

Talk-deck helper: `plot_eig_dist_vertical.py` is a tall 3-row λ overlay sized
for Keynote (native figsize, no `bbox_inches=tight`).

## Transfer / capacity gates (experimental)

Roadmap Track 1/2 diagnostics. They guide feature and selection work; they are
not production inference gates and do not replace DESI closure figures.

### G1 — GNN vs GBM on identical features

```bash
python workflows/sbi_inference/gate_g1_gnn_vs_gbm.py \
  --cache "${CACHE}" \
  --self-eval-npz /pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/abacus_self_linear_si/abacus_self_flowjax_preds.npz \
  --lambda-th 0.2
```

Trains `HistGradientBoostingRegressor` on the SI cache train split and scores
the same test rows as the Abacus self-eval npz. Adopt capacity / message-passing
work if `R²_GNN(λ₁) − R²_GBM(λ₁) > ~0.03`; otherwise the hand-crafted feature
set is likely information-limited. Asserts test-row and truth alignment between
cache and npz.

### G1.5 / G2 — RSD penalty and luminosity weighting

```bash
python workflows/sbi_inference/gate_g15_g2_rsd_luminosity.py \
  --ra 120 160 --dec 14.5 30.6 --zr 0.2 0.3 \
  --target-n 120000 --apertures-hmpc 3 7 10 14 \
  --lambda-th 0.2 --mass-log-th 13.0 --folds 5 --seed 42
```

Builds aperture features on one downsampled BGS-like Abacus cutsky sample in a
`{z-space, real-space} × {count, count+luminosity}` grid. Reports:

- **G1.5** RSD penalty ≈ `R²(real, count) − R²(z, count)` (upper bound on
  LOS-aware gains for this feature family);
- **G2** luminosity gain ≈ `R²(z, count+lum) − R²(z, count)` (GO if `>~0.03`).

Master cutsky path is hard-coded; needs `Z` and `Z_COSMO`, `R_MAG_ABS`, and
rs7 `LAMBDA1` truth. CPU only.

### A2 — mock vs DESI n(z)

```bash
python workflows/sbi_inference/measure_nz_mock_vs_desi.py \
  --ra 120 160 --dec 14.5 30.6 \
  --zmin 0.05 --zmax 0.55 --dz 0.01 \
  --out-dir /pscratch/sd/d/dkololgi/abacus/nz_comparison_YYYYMMDD
```

Same RA/Dec box for path1 mock parent and DESI BGS bright catalogue so shell
counts share solid angle. Explicitly counts and excludes mock sentinel phantoms
near `z ≈ 0.59` before writing the per-shell table, figure, and JSON.

## Figure publishing

During active runs, keep run directories as the source of truth. Mirror final
figures into the canonical browsable figure root with:

```bash
scripts/sync_figures_to_canonical.sh "${RUN_DIR}" desi_wedge_flowjax_linear_si
```

The destination defaults to
`$GRAPHWEB_CANONICAL_FIGURE_DIR/desi_wedge_flowjax_linear_si`, or
`/pscratch/sd/d/dkololgi/graphweb_desi/figures/desi_wedge_flowjax_linear_si`
when the environment variable is unset.

## Common pitfalls

- Do not mix the regression wedge and path1 fiberassign wedge for
  `--abacus-gnn-arrays`; array shapes can match while scaler constants do not.
- Do not run on legacy Mpc/h DESI wedge artifacts; inference intentionally aborts
  on those metadata.
- Do not use `--scale-invariant-features` with a baseline cache/model, or a
  baseline transform with a scale-invariant cache/model.
- `run_infer_desi_wedge_flowjax.sh` launches the baseline linear run; SI
  production needs an edited command with matching SI artifacts.
- `SFR <= 0` is treated as quenched in closure fractions and excluded from median
  sSFR panels.
- CIGALE and FastSpecFit property tables must not be mixed in one figure set
  without noting `sfr_source`; repoint `WEDGE_PARQUET` explicitly.
- Gate scripts (G1 / G1.5 / G2 / n(z)) hard-code Abacus/DESI paths and are
  Abacus-domain diagnostics — do not treat their GO/NO-GO prints as DESI VAC
  acceptance.
- `scripts/sync_figures_to_canonical.sh` syncs only top-level matching image /
  video / HTML files from the source directory (`rsync --exclude='*'`); nested
  artifacts need a separate sync command.
