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

Graph-construction constraints that this path inherits (row-order, no
north–south edges, public vs collaboration zall, cuGraph `--out-prefix`) are in
`RUNBOOK.md` under Catalog assembly and Gudhi/cuGraph.

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
  embeddings, and a posterior-sample subset (`--save-sample-subset`, default
  20,000 galaxies; full `[N,128,3]` samples are not written).
- `summary.json`: provenance, scaler statistics, class fractions, ordering
  checks, and Abacus reference fractions.

Post-hoc λ sort (default on): after converting flow samples to physical
eigenvalues, the script prints `true_ordering_violation_rate` (fraction of
samples that are not already `λ1≤λ2≤λ3`), then sorts each sample ascending
unless you pass `--no-sort`. Physical T-web eigenvalues are ordered; the
T-web *count* is order-independent, so the sort is a monotonic eval-time
correction, not a retraining. `summary.json` records `post_hoc_sort_applied`.
`--chunk-size` (default 512) only controls posterior sampling batches.

`summary.json` `reference_fractions.regression_desi` is a **frozen** Jraph
comparison (`void=0.246`, `wall=0.462`, `filament=0.266`, `cluster=0.026`),
not a live recompute. `plot_desi_wedge_flowjax.py --no-regression` drops those
bars. Do not treat them as the current Jraph product.

Unlike Jraph inference, this script **aborts** if DESI metadata
`coordinate_units` is `mpc/h`, `mpc_per_h`, or unknown, and it **does not**
force `JAX_PLATFORMS=cpu` (the checked-in launcher is a 1×A100 `salloc`).
Baseline (non-SI) runs also assert path1 fiberassign edge-scaler constants;
passing the Jraph regression-wedge NPZ fails that check even when shapes
match. GNN hparams come from the FlowJAX model pickle, not CLI.

### Domain-adaptation and production flags

The inference script exposes three transfer-diagnostic modes:

- `--edge-domain-adapt`: re-standardizes DESI scaled edge length and density
  contrast to training-like `N(0,1)`.
- `--node-domain-adapt`: re-standardizes DESI box-cox node features to
  training-like `N(0,1)`.
- `--scale-invariant-features`: uses per-graph-median normalized node and edge
  scale features. Use this only with a matching scale-invariant training cache
  and model.

SI node transform (must match the SI cache): divide columns
`[0, 2, 3, 4, 5, 6]` (`Degree`, `Density`, `Neigh Density`, `I_eig1`,
`I_eig2`, `I_eig3`) by their per-graph medians, then apply the cache
`PowerTransformer`. **Clustering (col 1) is left alone** — it is already a
dimensionless 0–1 fraction. SI edges divide `edge_length` (col 0) by the
per-graph median in `shared/abacus_gnn_parity.py` *before* log+z-score; the
path1 scaler-constant assert is skipped in SI mode.

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

There is no `--scale-invariant-features` flag: the Abacus graph is already
preprocessed inside the calibration cache. Pass the **matching** SI (or
baseline) model *and* cache together. The argparse default `--output-dir` is
`.../abacus_self_linear` (baseline); SI comparison plots need an explicit SI
directory as in the example. Samples are always sorted ascending (no
`--no-sort`). The script evaluates the cache **test** mask only.

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
- `--abacus-self-npz`: `abacus_self_flowjax_preds.npz` for the 3-way λ overlay
  (`eig_truth`, `lambda_mean`) and embedding PCA (`embeddings`). This file is
  **not** the class-fraction "Abacus NPE" series.
- `--abacus-classprob-npz`: class-fraction "Abacus NPE" bars. The NPZ must
  have 1-d keys `void`, `wall`, `filament`, `cluster` (probability arrays;
  the plot takes argmax). `abacus_self_flowjax_preds.npz` stores stacked
  `classprob` `[N,4]` and `hard_class` instead — passing it here KeyErrors
  on `void`. Without this flag the Abacus NPE bars are NaN even when
  `--abacus-self-npz` is set.
- Sky-box flags (`--ra-min` … `--z-max`) default to the expanded wedge
  (RA 120–160, Dec 14.5–30.6, z 0.2–0.3) and are used for the survey-edge
  panel. `--max-3d-points` / `--umap-points` cap optional HTML figures.

Common products include:

- `class_fractions_comparison.png`;
- `class_sky_map.png`;
- `posterior_width_sky_map.png`;
- `eigenvalue_distributions_3way.png` (the on-disk three-way λ overlay; some
  comments/print labels still say `domain_shift_overlay.png`);
- `embedding_pca.png` and optional `embedding_pca_3d.html`;
- optional UMAP and 3D HTML figures when dependencies are available.

For the cluster-deficit transfer story, see:

- `CLUSTER_DEFICIT_FIX_PLAN.md` (decision tree: Phase 0 domain correction →
  Route A scale-invariant features → optional mock densification);
- `plot_edge_scale.py` / `plot_cluster_recovery_bars.py` (graph-scale offset and
  recovery summary);
- `plot_eigenvalue_corner.py` / `plot_eig_dist_buildup.py` (λ joint structure).

### Cluster-deficit falsification index

These scripts recorded the June 2026 falsification of FoG / under-density /
training-extrapolation / dense-structure-shape as the *sole* cause of the DESI
cluster deficit. Most are CPU-only, take few or no CLI flags, and still
hard-code the **baseline** `desi_wedge_flowjax_linear` (or matching Abacus
self-eval) paths — retarget before comparing against SI production runs.

| Script | Intent |
| --- | --- |
| `plot_fog_los_alignment.py` | FoG test: major inertia eigenvector vs LOS (`\|cos θ\|`) on DESI vs Abacus. |
| `plot_training_coverage.py` | Baseline box-cox support: DESI vs Abacus train p99.9 / max per node feature. |
| `desi_abacus_coverage_report.py` | Same idea in **SI** model-input space (per-graph-median cols + cache PowerTransformer). CLI: `--si-cache`, `--desi-gnn-arrays`, `--preds-npz`. |
| `plot_shape_misclassification.py` | Dense-galaxy cluster rate vs local anisotropy; Oaxaca-style gap split. |
| `mmd_misspecification_check.py` | Unbiased RBF MMD² on GNN embeddings vs Abacus split-half floor. CLI: `--abacus-npz`, `--desi-npz`. |
| `plot_lambda_th_sweep.py` | Class fractions vs `λ_th` for Abacus truth / Abacus NPE / DESI NPE (+ Plotly morph). |
| Talk visuals | `build_skewer_animation.py`, `build_skewer_idealised.py`, `render_skewer_video.py`, `render_class3d_video.py`, `plot_eig_dist_buildup.py`, `plot_eig_dist_vertical.py` — see below. |

Related longer investigations (also non-canonical). Use them when reproducing a
specific SCIENCE_LOG note, not as weekly launch defaults.

| Script | Intent / constraint |
| --- | --- |
| `velocity_dispersion_precheck.py` | Phase-0 FoGAniso from Delaunay-neighbour LOS scatter. CLI: `--desi-preds-npz`, `--desi-gnn-arrays`, mock xyz/GNN/targets. CPU. |
| `velocity_dispersion_aperture_precheck.py` | Same kinematics at fixed apertures (~7 Mpc/h T-web scale). Neighbour Delaunay scale had no cluster signal; this tests the matched smoothing scale. |
| `velocity_dispersion_eigenvalue_precheck.py` / `threshold_limit_investigation.py` | Boundary-zone λ₁ predictability / metric artifact vs noise floor. Mock truth. |
| `smoothing_scale_investigation.py` / `plot_smoothing_scale_study.py` / `cluster_recovery_vs_smoothing.py` | Features fixed, T-web **target** smoothing varies (rs 6–24 Mpc/h cutsky FITS, hard-coded). High-memory CPU. |
| `mass_anchored_cluster_test.py` / `plot_mass_anchored_recovery.py` | Anchor “cluster” to `HALO_MASS` (column is **1e10 Msun/h**; `log10(M/[Msun/h]) = log10(HALO_MASS)+10`). Master cutsky is row-aligned with the rs catalogs. |
| `property_ceiling_ablation.py` | Redundancy of FastSpecFit properties vs 80-d GNN embedding (CLI: `--preds-npz`, `--closure-parquet`). DESI has no true env label, so this is the honest headroom proxy. |
| `probe_halo_mass_join.py` / `join_validate.py` | Abacus CompaSO `(FILE_NUM, BOX_INDEX)` indexing probes. **Not** the DESI FastSpecFit property join. Hard-coded path1 / halo_info paths. |

## Talk visuals (Keynote)

These scripts are presentation renderers, not inference. They consume
**FlowJAX** products (`desi_wedge_flowjax_preds.npz` keys `ra`/`dec`/`z`,
`hard_class`, `lambda_mean`, `lambda_std`, `classprob`). They KeyError on
Jraph `desi_wedge_index_and_preds.npz` (`cls`, no `hard_class`). The inverse
trap is `workflows/visualization/visualize_desi_wedge_cweb_3d.ipynb`, which
is Jraph-only.

### Real DESI skewer

`build_skewer_animation.py` draws a radial pencil beam through the richest
inferred cluster:

- Sightline: galaxy with max `classprob[:, 3]`, then a `classprob`-weighted
  RA/Dec centroid of neighbours within `--anchor-r` (default 15 Mpc).
- Beam: angular radius `--theta-deg` (default 0.6).
- Per-frame λ densities: mixture of per-galaxy `lambda_mean` / `lambda_std`
  along the beam. The module docstring still mentions `lambda_samples_subset`;
  the code does not read that key.
- Default `--out`: `{preds parent}/skewer_posterior_animation_real.html`.

```bash
python workflows/sbi_inference/build_skewer_animation.py \
  --preds-npz "${RUN_DIR}/desi_wedge_flowjax_preds.npz" \
  --theta-deg 0.6 --anchor-r 15
```

### Idealised skewer (no catalogue)

`build_skewer_idealised.py` synthesises a Gaussian density bump so λ3, λ2, λ1
cross `λ_th=0.2` at d≈0.33, 0.60, 0.86 (void→wall→filament→cluster). It
imports the HTML template from `build_skewer_animation.py` and takes **no**
preds NPZ. Default `--out` writes
`.../desi_wedge_flowjax_linear_si/skewer_idealised.html` even if that SI run
directory does not exist yet.

### HTML → mp4 / gif

`render_skewer_video.py` extracts the embedded payload (`var D={...};var F=`)
and re-renders with matplotlib. Default `--out` is `.mp4` (needs ffmpeg on
PATH); pass `--out ....gif` for Pillow. Works for both the real and idealised
HTML files.

```bash
python workflows/sbi_inference/render_skewer_video.py \
  "${RUN_DIR}/skewer_posterior_animation_real.html"
```

### Rotating class-3D GIF

`render_class3d_video.py` rebuilds the Plotly class-3D cloud as a matplotlib
fly-through GIF (void→wall→filament→cluster reveal, then rotation). Default
`--preds` is the SI FlowJAX npz. It keeps **all** clusters and subsamples the
other classes (`--n-points`, default 28,000). Default output:
`{preds parent}/class_3d_flythrough.gif`.

### Talk λ overlays

- `plot_eig_dist_buildup.py`: two PNGs with **shared bins** for successive
  slides (`eig_dist_buildup_abacus_only.png`, then
  `eig_dist_buildup_with_desi.png`).
- `plot_eig_dist_vertical.py`: tall 3-row PDF/PNG, native figsize 7.0×7.3 in,
  **no** `bbox_inches=tight` so Keynote point sizes match the saved file.
  Writes both `eig_dist_vertical_abacus_only` and `eig_dist_vertical_3way`.
  Defaults to SI Abacus-self / DESI-preds / figure paths.

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

The join recovers `TARGETID` through `global_node_id` as a **row index** into
the source maglim FITS (same order as the Gudhi graph). Default `--preds` is
the SI production npz; override it if you are not on that run. Default
`--fastspec-glob` is `fastspec-loa-main-bright-nside1-hp*.fits` under
`DESI_FASTSPEC_CATALOGS_DIR`. Hemisphere-split graphs duplicate some
`TARGETID`s (N/S copies with different neighbourhoods); the join averages
posterior products, sets `is_dup`, and recomputes `hard_class` from the
averaged class probabilities.

Rest-frame colour in the parquet is FastSpecFit `ABSMAG01_SDSS_G − ABSMAG01_SDSS_R`
(z=0.1 SDSS). That is **not** the LSS k+e `ABSMAG_RP1` from
`desi_absmag_kcorr.py`. Do not mix the two magnitude systems in one figure.

Writes `desi_wedge_env_props.parquet` plus `join_report.json`. The report
includes FastSpecFit match fraction, duplicate-target count, valid-sSFR
fraction, and an RA/Dec cross-check.

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
  FlowJAX baseline mode asserts path1 constants and aborts; Jraph inference
  silently refits whatever NPZ you pass.
- Do not run FlowJAX on legacy Mpc/h DESI wedge artifacts; it aborts. Jraph
  inference only **warns** and continues — a silent unit mix-up is possible.
- Do not use `--scale-invariant-features` with a baseline cache/model, or a
  baseline transform with a scale-invariant cache/model. SI node medians skip
  Clustering (col 1).
- `run_infer_desi_wedge_flowjax.sh` launches the baseline linear run; SI
  production needs an edited command with matching SI artifacts.
- `infer_abacus_self_flowjax.py` default output dir is `abacus_self_linear`;
  SI three-way plots need an explicit SI directory and the matching SI cache.
- `--no-sort` keeps raw flow λ order; default DESI inference and Abacus self
  both sort samples ascending after recording the true violation rate.
- `summary.json` `regression_desi` fractions are frozen Jraph numbers, not a
  live Jraph re-run. Use `--no-regression` unless you want those bars.
- `SFR <= 0` is treated as quenched in closure fractions and excluded from median
  sSFR panels.
- CIGALE and FastSpecFit property tables must not be mixed in one figure set
  without noting `sfr_source`; repoint `WEDGE_PARQUET` explicitly.
- Gate scripts (G1 / G1.5 / G2 / n(z)) hard-code Abacus/DESI paths and are
  Abacus-domain diagnostics — do not treat their GO/NO-GO prints as DESI VAC
  acceptance.
- Cluster-deficit falsification helpers often hard-code baseline
  `desi_wedge_flowjax_linear` products; retarget before SI comparisons.
- `scripts/sync_figures_to_canonical.sh` syncs only top-level matching image /
  video / HTML files from the source directory (`rsync --exclude='*'`); nested
  artifacts need a separate sync command.
- `join_validate.py` is an Abacus CompaSO index check, not the DESI property
  join (`build_desi_wedge_property_join.py`).
- `HALO_MASS` on the master cutsky is in **1e10 Msun/h**; mass-anchored scripts
  add +10 after `log10`.
- Do not point `GalaxyCatalog` or `load_catalog.py` products at this path; SBI
  needs the maglim FITS with `TARGET_RA`/`TARGET_DEC` in the same row order as
  the Gudhi graph.
- Property-join `gr` is FastSpecFit `ABSMAG01_SDSS_*`, not LSS `ABSMAG_RP1`.
- The 3D notebook `visualize_desi_wedge_cweb_3d.ipynb` reads Jraph
  `desi_wedge_index_and_preds.npz`, not `desi_wedge_flowjax_preds.npz`.
- Skewer / class-3D GIF scripts are the inverse: they need FlowJAX
  `hard_class` / `lambda_mean` / `lambda_std` / `classprob`. Do not point
  them at the Jraph npz.
- `--abacus-self-npz` does not fill class-fraction "Abacus NPE" bars.
  `--abacus-classprob-npz` needs per-class 1-d keys, not stacked `classprob`.
- `build_skewer_animation.py` mixes `lambda_mean`/`lambda_std`; it does not
  read `lambda_samples_subset` despite the module docstring.
- `render_skewer_video.py` default mp4 needs ffmpeg; a missing encoder fails
  after the HTML already exists. `build_skewer_idealised.py` default `--out`
  points at the SI run directory even when that run has not been created.
