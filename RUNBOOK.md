# GraphWeb_DESI Runbook

This runbook lists validated workflow entrypoints and launch commands.

## Environment setup

Recommended shell setup for GraphWeb workflows:

```bash
source ~/.bashrc
conda activate cosmic_env
desienv
```

If `desimodel` import remains unavailable in your shell context, use the
canonical script paths and environment used in validation notes. The catalog
workflow also includes a hardcoded `desimodel` fallback path.

The cuGraph feature-export step needs the `rapids-gnn` GPU env; GAT, Gudhi
graph build, Jraph, and FlowJAX inference use `cosmic_env` (see
`~/.claude/CLAUDE.md`). The absolute-magnitude diagnostic below uses the DESI
software stack instead.

## Canonical workflows

### Catalog assembly

There are two DESI catalogs. Do not mix them.

#### Low-z FastSpecFit product (GAT stack only)

Canonical:

```bash
python workflows/catalog/load_catalog.py --help
```

Compatibility wrapper:

```bash
python load_catalog.py --help
```

Example run (`z ≈ 0.01–0.06` only):

```bash
python workflows/catalog/load_catalog.py \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/catalogs
```

`--out-dir` defaults to `GRAPHWEB_CATALOG_DIR` (pscratch). Do not write these
FITS products under the repo or `$HOME` — home previously hit its 40 GiB quota.
Resolve paths via `shared/config_paths.py`:

- `GRAPHWEB_CATALOG_DIR`
- `GRAPHWEB_CATALOG_PATH` (`loa-combined-lowz.fits`)
- `GRAPHWEB_CATALOG_ZFLAGS_PATH`
- `GRAPHWEB_CATALOG_FASTSPEC_PATH`

Selection in code: `SPECTYPE==GALAXY` and `--z-min/--z-max` (defaults 0.01–0.06),
then a left join to the redshift catalog keeping `ZCAT_PRIMARY==True`, then
optional Legacy Survey Tractor photometry. The FastSpecFit file list is a
hard-coded `FASTSPEC_CATALOGS` array in `load_catalog.py` (not a directory
glob). `--zcat-file` defaults to `DESI_ZCAT_FILE` in `shared/config_paths.py`
(`/global/cfs/cdirs/desi/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits`).
`desimodel.footprint.radec2pix` is required for the Tractor join; the script
tries a Perlmutter desiconda fallback if `desienv` is not active.

This low-z catalog does **not** overlap the GraphWeb-BGS VAC redshift range
(`0.15–0.55`) and cannot supply `ABSMAG` for that work. Columns include `RA` /
`DEC` (FastSpecFit METADATA), which is what `GalaxyCatalog` expects.

#### Bright BGS maglim product (Gudhi / Jraph / SBI)

```bash
python workflows/catalog/build_bgs_maglim_catalog.py \
  --out-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits
```

Streaming filter of the LOA `zall` `ZCATALOG` HDU (2e6-row chunks; no full
in-memory load). Kept rows:

- `ZWARN == 0`
- `DELTACHI2 >= 25`
- `SPECTYPE == GALAXY`
- `BGS_TARGET` has a BGS_BRIGHT bit (`BGS_BRIGHT` / `_NORTH` / `_SOUTH`) unless
  `--no-bright-only` (then any nonzero `BGS_TARGET`)

No redshift cut and no stellar-mass cut. Output columns are the minimal zall
set (`TARGETID`, `SURVEY`, `PROGRAM`, `TARGET_RA`, `TARGET_DEC`, `Z`, `ZWARN`,
`DELTACHI2`, `SPECTYPE`, `BGS_TARGET`). Existing `--out-path` refuses to
overwrite unless `--overwrite`.

Default `--zall-path` is the **public** DR2 copy
`/global/cfs/cdirs/desi/public/dr2/spectro/redux/loa/zcatalog/v1/zall-pix-loa.fits`.
That is **not** `DESI_ZCAT_FILE` used by `load_catalog.py`. Pass `--zall-path`
explicitly if you need the collaboration redux file.

### Graph inference pipeline

Canonical:

```bash
python workflows/graph_inference/graph_catalog.py --help
```

Compatibility wrapper:

```bash
python graph_catalog.py --help
```

Example no-plot cache-only probe:

```bash
export GRAPHWEB_CACHE_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/cache
export GRAPHWEB_VAC_OUTPUT_PATH=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/DESI_BGS_PRERELEASE_VAC.pkl
python workflows/graph_inference/graph_catalog.py --cache-mode cache-only --no-summary-plot
```

Defaults still write under the **repo** (`GRAPHWEB_CACHE_DIR={repo}/cache`,
`GRAPHWEB_VAC_OUTPUT_PATH={repo}/DESI_BGS_PRERELEASE_VAC.pkl`).
`GRAPHWEB_CANONICAL_CACHE_DIR` on pscratch is unused by this script — set the
env vars (or `--cache-dir` / `--vac-output-path`) before a `rebuild`, or home
quota fills the same way the low-z FITS catalogs did.

`--no-summary-plot` skips the histogram but still constructs a TNG300 `cat`
from `--reference-catalog-path` (`TNG_REFERENCE_CATALOG_PATH`). A missing TNG
tree fails even a cache-only no-plot probe.

`--first-moment-matching` applies the Illustris training `features_scaler.pkl`
instead of fitting a Box-Cox `PowerTransformer` on DESI, then **still**
subtracts DESI column means. The `--help` fast-exit usage line omits this flag;
it exists on the real argparse parser. Leave it off unless you are deliberately
matching the TNG scaler.

Rebuild constructs an Illustris `network(masscut=9.0, from_DESI=True)` and
feeds `SimpleGAT(input_dim=10, output_dim=4, num_heads=4)`. Those 10 GAT node
features are **not** the 7 Abacus-style columns used by Jraph/SBI
(`Degree, Clustering, Density, Neigh Density, I_eig1, I_eig2, I_eig3`). The
GAT alpha/Delaunay NetworkX cache is a different graph from the Gudhi/cuGraph
Mpc wedge.

VAC pickle columns written onto the zcat snapshot: `GAT_ENV` (argmax class
0–3) and `GAT_VOID_PROB` / `GAT_WALL_PROB` / `GAT_FILAMENT_PROB` /
`GAT_CLUSTER_PROB`. Cache filenames are not symmetric: alpha uses
`DESI_NETWORKalpha_zcat.pt`; Delaunay uses `DESI_NETWORK_delaunay_zcat.pt`.
`--vac-output-path` must include a directory component (`write_vac` calls
`os.makedirs(dirname(...))`). The histogram uses `scienceplots` colours, not
`shared.plot_style.COSMIC_WEB_COLORS`.

If model checkpoint path differs from config default:

```bash
python workflows/graph_inference/graph_catalog.py --model-path /path/to/trained_gat_model.pth
```

### Gudhi/cuGraph/Jraph wedge inference

This path is for DESI wedge eigenvalue regression with an Abacus-trained Jraph
model. It is not a replacement command for the all-sky GAT VAC workflow above.
Run it on Perlmutter from an environment with DESI, Gudhi/cuGraph where needed,
JAX/Haiku/Jraph, and access to the Illustris training repo.

Build the bright BGS catalog if it does not already exist (selection and
`--zall-path` defaults are under Catalog assembly above):

```bash
python workflows/catalog/build_bgs_maglim_catalog.py \
  --out-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits
```

Build the full Mpc-parity graph and export graph features:

```bash
python workflows/graph_construction/build_desi_bgs_gudhi_graph.py \
  --catalog-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc \
  --alpha-sq inf \
  --coord-units mpc

python workflows/graph_construction/desi_graph_features_cugraph.py \
  --metadata-path /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_delaunay_metadata.json \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc \
  --out-prefix desi_bgs_cugraph
```

Construction constraints (from the scripts, not the launchers):

- The maglim FITS must stay in the same row order from catalog → Gudhi graph →
  cuGraph arrays → wedge subset. The wedge cutter maps `TARGET_RA` /
  `TARGET_DEC` / `Z` by catalog row index.
- Gudhi splits Galactic hemispheres (`b>0` vs `b<=0`) and builds one
  AlphaComplex each. There are **no north–south edges**. `--alpha-sq inf`
  (default) is full Delaunay; a finite value writes `desi_alpha_*` instead of
  `desi_delaunay_*`.
- `--coord-units mpc` (default) is Abacus training parity. `mpc_per_h` is the
  deprecated ×h convention.
- Full-hemisphere AlphaComplex on ~10⁷ galaxies is expensive. Use
  `--max-points-per-hemi N` only for smoke tests (random subsample, `--seed`).
  The wedge cutter's `--max-nodes N` is a different cap: it keeps the **first
  N selected catalog rows**, not a random subsample.
- cuGraph feature export must run in `rapids-gnn` on a GPU node. The script
  default `--out-prefix` is `desi_delaunay_cugraph`; the validated expanded
  wedge uses `--out-prefix desi_bgs_cugraph` to match
  `JRAPH_INPUTS_expanded_wedge.txt`. `--skip-clustering` / `--skip-inertia`
  zero those columns (not production). `--write-parquet` is huge on the full
  graph. `--debug-inertia-parity` compares two inertia implementations on a
  random subgraph and **exits**; `--max-edges-export N` truncates the exported
  edge table.
- On Slurm, `conda activate` inside `srun bash -lc` is not enough: unset
  `PYTHONPATH`/`PYTHONHOME`, set `PYTHONNOUSERSITE=1`, and call the env’s
  absolute `python` (see `.cursor/rules/conda-env-srun-python-path.mdc`). A
  login-node `import gudhi` does not prove the compute node will find it.

Subset the expanded Abacus-parity DESI wedge and run inference:

```bash
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

python workflows/graph_construction/subset_desi_graph_wedge.py \
  --graph-metadata /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_delaunay_metadata.json \
  --catalog-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --parent-gnn-arrays /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_bgs_cugraph_gnn_arrays.npz \
  --parent-gnn-metadata /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_bgs_cugraph_gnn_metadata.json \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph \
  --out-prefix desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc \
  --ra-min 120 --ra-max 160 --dec-min 14.5 --dec-max 30.6 --z-min 0.2 --z-max 0.3

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

`INFER_DIR` in `workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt` is the full
intended run directory, so the command passes its parent as `--output-dir` and
the basename as `--run-name`. Jraph inference uses `ILLUSTRIS_ROOT` to find
`shared/graph_net_models.py`; the GAT classifier uses `ILLUSTRIS_REPO_ROOT` from
`shared/config_paths.py`.

Validated one-liner after the expanded Mpc wedge products exist:

```bash
bash workflows/catalog/run_infer_expanded_wedge_mpc.sh
```

That wrapper sources `JRAPH_INPUTS_expanded_wedge.txt`, forces `cosmic_env`
Python with a clean `PYTHONPATH`, sets `JAX_PLATFORMS=cpu`, and writes under
`/pscratch/sd/d/dkololgi/graphweb_desi/inference_outputs/${INFER_RUN_NAME}`.
The Python script itself also defaults `JAX_PLATFORMS=cpu` and clears
`CUDA_VISIBLE_DEVICES` unless you already set a platform — unlike FlowJAX,
this path will not use a GPU unless you override that.

The inference script loads `shared/abacus_gnn_parity.py` from **this** repo
(bidirectional edges, `1/density_contrast` on the reverse, log+z-score of
length and density-contrast) and `shared/graph_net_models.py` from
`ILLUSTRIS_ROOT`. Do not put GraphWeb_DESI on `sys.path` ahead of Illustris or
the wrong `shared` package wins. The script loads that parity helper by file
path (`importlib`) so Illustris `shared/` does not shadow it.

Checkpoint and architecture (from the code, not the model pickle):

- `--abacus-run-dir` must contain `checkpoints/best_checkpoint.json` with a
  `path` field; that pickle is the weights file.
- `--latent-size/--num-heads/--num-passes/--dropout` default to `96/8/8/0.15`
  (manifest `LATENT_SIZE` etc.). They are **not** read from the checkpoint.
  Mismatch with training is a silent architecture error.
- `output_dim` is hardcoded to **15**. This companion is the 15-d halo_xcom
  regression path, not a 3-d increment NPE.

Units and scaler guards differ from FlowJAX:

- Jraph **warns** if DESI `coordinate_units` is `mpc/h`, `mpc_per_h`, or
  unknown, then continues. FlowJAX **aborts**. Do not rely on Jraph to catch
  a leftover Mpc/h wedge.
- Jraph **fits** the edge scaler from whatever `--abacus-gnn-arrays` you pass.
  FlowJAX baseline mode asserts path1 fiberassign constants. Mixing the
  regression-wedge NPZ into Jraph is a silent scaler bug.

Validated 15-d caches store **ordered linear increments** after
`target_scaler.inverse_transform`: `v1=λ1`, `v2=λ2−λ1`, `v3=λ3−λ2`. Reconstruct
with `λ2=λ1+max(v2,ε)`, `λ3=λ2+max(v3,ε)`. Do not treat `v2`/`v3` as λ₂/λ₃.
If the calibration cache has `node_feature_scaler` (training
`--power-scale-node-features`), DESI raw `x` must be
`scaler.transform(x+1e-6)` before the forward pass; skipping that step
collapses class fractions (~100% void). The script warns and continues if the
scaler key is missing.

Outputs under `--output-dir/--run-name/` include `preds_scaled_15d.npy`,
`preds_raw_15d.npy`, `preds_lambda123.npy`, `class_lambda_thr0.2.npy`,
`summary.json`, histogram/bar PNGs, and `desi_wedge_index_and_preds.npz`
(`global_node_id`, `ra/dec/z`, `lambda1/2/3`, `cls`). The visualization
notebook and `JRAPH_INPUTS_expanded_wedge.txt` `INFER_NPZ` point at that npz,
not the FlowJAX `desi_wedge_flowjax_preds.npz`.

Do **not** use these as substitutes for the Gudhi NPZ path:

- `workflows/jraph_inference/build_desi_wedge_jraph_cache.py` — induced wedge
  from the GAT Delaunay NetworkX cache under `{repo}/cache/`, not cuGraph
  arrays. Reverse density-contrast is stored as `-dc`, which does not match
  Abacus `1/dc` in `shared/abacus_gnn_parity.py`.
- `workflows/jraph_inference/experiment_desi_feature_parity.py` — frozen
  2026-05-29 one-off against **Mpc/h** DESI arrays
  (`..._bright_from_fullgraph/`, no `_mpc` suffix). It diagnosed the unit
  mismatch; do not re-run it against current Mpc products.

Path-1 Jraph 3D comparison notebook:
`workflows/visualization/visualize_desi_wedge_cweb_3d.ipynb`. Edit **cell 1**
only (paths). It loads Jraph `INFER_DIR/desi_wedge_index_and_preds.npz` plus
the expanded Mpc DESI GNN arrays and Abacus training-wedge truth. Interactive
HTML lands under `INFER_DIR`; download it off Perlmutter rather than saving
huge Plotly outputs back into the notebook. This is not the FlowJAX plotter.
`workflows/visualization/path1_desi_wedge_inference_summary.md` is a NERSC-only
symlink into pscratch and is dangling off Perlmutter.

#### Perlmutter Slurm chain (graph → features → wedge)

`workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh` is the compute-node
launcher for the Gudhi / cuGraph / wedge stages (never run those on login).
Stages:

| `STEP` | Role | Env / hardware |
| --- | --- | --- |
| `2` (default) | Full Gudhi graph (`build_desi_bgs_gudhi_graph.py`) | `cosmic_env`, CPU, long walltime |
| `2b` | cuGraph GNN feature export | `rapids-gnn`, GPU |
| `3` | RA/Dec/z wedge subset | `cosmic_env`, CPU |

Chain from login (lightweight submit only):

```bash
SUBMIT_CHAIN=1 bash workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh
```

Or submit stages explicitly with dependencies (see the script header for the
full `sbatch` flags). Logs land under
`/pscratch/sd/d/dkololgi/graphweb_desi/logs/desi_bgs_bright_pipeline/`.

**Important path caveat:** the checked-in launcher still hard-codes the older
narrow bright-wedge directory names and sky cut
(`gudhi_hemi_full_alphasq_inf_seed42_bright`, RA 120–140 / Dec 16.5–26.7 /
z 0.25–0.3, no `_mpc` suffix). Graph build uses the script default
`--coord-units mpc`, but those output paths are **not** the validated expanded
Mpc products in `JRAPH_INPUTS_expanded_wedge.txt`. For Abacus-parity Jraph/SBI
work, prefer the expanded-wedge commands in this runbook (or edit the launcher
paths/`--ra-*`/`--dec-*`/`--z-*` before submitting). Step `2b` still requires
the `rapids-gnn` GPU env regardless of which wedge geometry you choose.

### FlowJAX/SBI DESI wedge posterior inference

This path runs an Abacus-trained FlowJAX neural posterior estimator on the DESI
wedge. It returns per-galaxy posterior eigenvalue summaries and class
probabilities, not just a hard class or point-regression output. Use it after the
Mpc-parity Gudhi/cuGraph wedge products above exist.

Detailed notes and caveats:

```bash
less workflows/sbi_inference/README.md
```

Validated Perlmutter launcher (baseline linear run — no SI flag):

```bash
bash workflows/sbi_inference/run_infer_desi_wedge_flowjax.sh
```

Manual baseline command shape:

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
`--scale-invariant-features`. SI node medians skip Clustering (col 1). For
transfer diagnostics only, the script also exposes `--edge-domain-adapt` and
`--node-domain-adapt`. Posterior samples are sorted ascending unless
`--no-sort`; `summary.json` `regression_desi` fractions are frozen Jraph
numbers. Depth: `workflows/sbi_inference/README.md`.

Generate DESI diagnostics after inference:

```bash
RUN_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear_si

python workflows/sbi_inference/plot_desi_wedge_flowjax.py \
  --preds-npz "${RUN_DIR}/desi_wedge_flowjax_preds.npz" \
  --summary-json "${RUN_DIR}/summary.json" \
  --output-dir "${RUN_DIR}" \
  --no-regression
```

`--abacus-self-npz` fills the 3-way λ overlay / embedding PCA from
`abacus_self_flowjax_preds.npz`. It does **not** fill the class-fraction
"Abacus NPE" bars; that series needs `--abacus-classprob-npz` with 1-d keys
`void`/`wall`/`filament`/`cluster`. Without it those bars are NaN.

Talk-visual / Keynote renderers (`build_skewer_animation.py`,
`build_skewer_idealised.py`, `render_skewer_video.py`,
`render_class3d_video.py`, `plot_eig_dist_vertical.py`) are FlowJAX-only and
are documented in `workflows/sbi_inference/README.md`. Do not point them at
Jraph `desi_wedge_index_and_preds.npz`.

Join FastSpecFit properties and make closure figures:

```bash
python workflows/sbi_inference/build_desi_wedge_property_join.py \
  --preds "${RUN_DIR}/desi_wedge_flowjax_preds.npz"

python workflows/sbi_inference/plot_property_environment_closure.py \
  --table "${RUN_DIR}/desi_wedge_env_props.parquet" \
  --outdir "${RUN_DIR}"
```

Mirror final top-level figures into the canonical figure root:

```bash
scripts/sync_figures_to_canonical.sh "${RUN_DIR}" desi_wedge_flowjax_linear_si
```

Optional property-science path (CIGALE-HZ mass/SFR swap + SFMS / M*–colour /
mass-controlled figures) and Abacus-domain transfer gates (G1, G1.5/G2, n(z))
are documented in `workflows/sbi_inference/README.md`. They are not part of the
canonical inference launcher.

### DESI absolute-magnitude support diagnostic (experimental)

`workflows/sbi_inference/desi_absmag_kcorr.py` evaluates whether a luminosity
feature is worth further parity work. It joins BGS photometry by `TARGETID`,
calls the official DESI LSS `add_dered_flux` and `add_ke` implementations, and
writes k+e-corrected `ABSMAG_RP1`. It then prints marginal magnitude and colour
statistics beside a fixed raw Abacus cut-sky sample.

Run this in the DESI environment, not `cosmic_env`:

```bash
source /global/common/software/desi/desi_environment.sh main

python workflows/sbi_inference/desi_absmag_kcorr.py \
  --wedge /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --full /global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full.dat.fits \
  --out /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_absmag_rp1_gate_sub.fits \
  --zlo 0.15 --zhi 0.55 \
  --max-rows 400000 --seed 42
```

The `--lsscode` checkout must contain both `LSS/py/LSS` and `DESI_ke`; its
default is the complete shared NERSC checkout used when this workflow was
developed. The input named `--wedge` must provide `TARGETID`, `Z`,
`TARGET_RA`, and `TARGET_DEC`; `--full` supplies `FLUX_G`, `FLUX_R`, and the
corresponding Milky Way transmissions. The current script does not apply an
RA/Dec cut to the DESI input despite the historical option name.

Operational constraints:

- The default 400,000-row sample is selected reproducibly **after** the positive
  flux and half-open `zlo <= Z < zhi` cuts, and **before** `add_ke`. That order
  is required so the support gate fits an interactive wall: `add_ke` is
  single-threaded and was projected ~4 h over the full ~8.7 M joined rows.
  `--max-rows 0` keeps every z-selected row.
- The output is a support-study subsample, not a deployment catalogue. A full
  catalogue `ABSMAG` product is a separate longer run, only needed if the gate
  looks promising.
- The FITS output is written before the Abacus comparison runs. A later
  comparison failure does not invalidate the k+e columns, but it does mean the
  diagnostic did not complete.
- The Abacus path and its fixed sky/magnitude selection are currently
  hard-coded. There is no CLI pass threshold or nonzero exit status for a
  support mismatch.
- On `setuptools>=81`, the script shims a minimal `pkg_resources` module before
  importing DESI_ke (`smith_kcorr` / `findfile` still import it even though the
  call site is dead). Do not patch the shared read-only LSS checkout for this.

Interpret the printed offsets as exploratory diagnostics only. The present
comparison uses different DESI and Abacus footprints and unmatched redshift
distributions, and it does not prove that DESI `ABSMAG_RP1` and Abacus
`R_MAG_ABS` share the same band, reference-redshift, distance-modulus, and
evolution-correction conventions. Before adding luminosity to production
inference, compare the exact training and inference selections in redshift
bins (or with matched `n(z)`), establish the magnitude convention explicitly,
and define quantitative acceptance thresholds.

### Utility workflows

`workflows/utilities/galaxy_catalog.py` is a GAT-stack helper, not a launcher
for the maglim/Gudhi path. `GalaxyCatalog(PATH, LOGMSTAR=9.)` filters
`ZWARN==0`, `DELTACHI2>=25`, `LOGMSTAR>=9`, `BGS_TARGET!=0`, converts `RA`/`DEC`
to Planck18 comoving Mpc, and reorders the table to `[north, south]`. Running
the script as `__main__` loads `GRAPHWEB_CATALOG_PATH` (`loa-combined-lowz.fits`).
It will KeyError on the maglim catalog (`TARGET_RA`/`TARGET_DEC`, no
`LOGMSTAR`). `--help` is a stub (`Usage: python galaxy_catalog.py`) and does
not list flags.

`workflows/utilities/investigate_edges.py` is GAT-cache QA: it loads
`{GRAPHWEB_CACHE_DIR}/DESI_alpha_geom.pt` and prints edge-length stats. It does
not read Gudhi/cuGraph arrays. `--help` is the same stub pattern.

Canonical:

```bash
python workflows/utilities/galaxy_catalog.py --help
python workflows/utilities/investigate_edges.py --help
```

Compatibility wrappers:

```bash
python galaxy_catalog.py --help
python investigate_edges.py --help
```

### P12-A Loa FMPE VAC

Provisional halo48 posterior on Loa BGS. Not a replacement for the GAT pickle
or the Jraph/FlowJAX wedge npz. Column legend, quality bits, and the script
map: `docs/p12a_vac_operations.md`.

Canary (file path required; `prepare` will not overwrite `INPUTS_READY.json`):

```bash
PY=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
ROOT=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1
$PY workflows/p12a_vac/loa_trial.py prepare  --root "$ROOT"
$PY workflows/p12a_vac/loa_trial.py finalize --root "$ROOT"
$PY workflows/p12a_vac/loa_trial.py infer    --root "$ROOT"
```

Full survey (`plan` will not overwrite `PLAN.json`; workers verify the frozen
copy under `$ROOT/source`, not the git checkout):

```bash
ROOT=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1
$PY workflows/p12a_vac/production.py plan   --root "$ROOT"
$PY workflows/p12a_vac/production.py worker --root "$ROOT" --part "$SLURM_ARRAY_TASK_ID"
$PY workflows/p12a_vac/production.py merge  --root "$ROOT"
```

`workflows/p12a_vac/run_production.sh` is the Perlmutter array/merge wrapper.
It sets `ILLUSTRIS_ROOT` to the frozen candidate snapshot. Completion without
`science_release_ready: true` is still not a science release. Alignment
scripts under `workflows/p12a_vac/alignment_*.py` are read-only diagnostics.

### NERSC helper scripts

These are login-node sanity checks, not pipeline stages:

- `scripts/where_am_i.sh` — host, Slurm env, HOME/PSCRATCH.
- `scripts/check_quota.sh` — `myquota` if present, else `df` on HOME/PSCRATCH.
- `scripts/desi_prods.sh` — lists top-level `/global/cfs/cdirs/desi` trees.
- `scripts/sync_figures_to_canonical.sh <src_dir> <run_name>` — top-level
  png/pdf/html/svg/gif/mp4 only (`rsync --exclude='*'`).
- `templates/job.sbatch` — minimal Perlmutter template (`srun` required).

## Troubleshooting

| Symptom | Likely cause / fix |
| --- | --- |
| Home disk full / cannot commit | Catalogs **and** GAT cache/VAC defaults belong on pscratch. `GRAPHWEB_CACHE_DIR` / `GRAPHWEB_VAC_OUTPUT_PATH` still default to the repo; `GRAPHWEB_CATALOG_DIR` already points at pscratch. |
| GAT `--no-summary-plot` still fails | `graph_catalog.py` always constructs a TNG300 `cat` from `TNG_REFERENCE_CATALOG_PATH` before inference. |
| Jraph ~100% void | Calibration cache missing / unused `node_feature_scaler` while training used `--power-scale-node-features`. |
| Jraph λ₂/λ₃ look like increments | 15-d caches store `v2=λ2−λ1`, `v3=λ3−λ2`; reconstruct before classifying. |
| Jraph `shared.*` import is GraphWeb_DESI | `ILLUSTRIS_ROOT` must precede this repo on `sys.path` (`run_infer_expanded_wedge_mpc.sh` already does). |
| GAT-cache Jraph pickle vs Abacus | `build_desi_wedge_jraph_cache.py` is not the Gudhi NPZ path; reverse density-contrast convention differs. |
| FlowJAX abort on coordinate units | Rebuild DESI wedge with `--coord-units mpc`; do not use legacy Mpc/h products. |
| SI vs baseline mismatch | Match `--scale-invariant-features` to the training cache/model; the checked-in launcher is baseline-only. |
| Jraph ran on an Mpc/h wedge anyway | Jraph only **warns** on `coordinate_units`; FlowJAX aborts. Check `summary.json` `desi_coordinate_units`. |
| Jraph on GPU still uses CPU | Script/wrapper default `JAX_PLATFORMS=cpu`. Override only if you intend a GPU run. |
| Jraph checkpoint loads but predictions look random | CLI hparams (`96/8/8/0.15`) must match training; they are not in the pickle. `output_dim=15` is hardcoded. |
| Edge scaler looks wrong | FlowJAX: pass path1 fiberassign arrays (asserted). Jraph: it will happily fit a regression-wedge NPZ — check the path. |
| Viz notebook KeyError / missing npz | Point cell 1 at Jraph `desi_wedge_index_and_preds.npz`, not FlowJAX `desi_wedge_flowjax_preds.npz`. |
| Skewer / class-3D GIF KeyError on `hard_class` | Those scripts need FlowJAX `desi_wedge_flowjax_preds.npz`. The 3D notebook is Jraph-only; this is the inverse trap. |
| Class-fraction "Abacus NPE" bars are NaN | `--abacus-self-npz` is the 3-way λ overlay. Pass `--abacus-classprob-npz` with 1-d keys `void`/`wall`/`filament`/`cluster`; the self-eval npz's stacked `classprob` is the wrong schema. |
| `render_skewer_video.py` fails after HTML exists | Default output is mp4 and needs ffmpeg. Use `--out ....gif` for Pillow. |
| Wedge smoke test looks spatially biased | `--max-nodes` keeps the first N catalog-order rows. Gudhi `--max-points-per-hemi` is the random subsample. |
| Property `g-r` disagrees with `ABSMAG_RP1` | Join uses FastSpecFit `ABSMAG01_SDSS_*` (z=0.1); kcorr writes LSS `ABSMAG_RP1`. |
| IBM Plex fallback / DejaVu plots | Fonts are gitignored under `assets/fonts/`. Set `GRAPHWEB_FONT_DIR` or run the PLOT_STYLE_GUIDE §2.3 curl. |
| `desi_absmag_kcorr.py` `pkg_resources` ImportError | Use the repo script (it shims `pkg_resources`); or pin setuptools `<81` only if you must call DESI_ke outside this wrapper. |
| `add_ke` never finishes in interactive QOS | Keep the default `--max-rows 400000` z-cut subsample; full-catalogue k+e is a separate long job. |
| Figures missing after sync | `sync_figures_to_canonical.sh` only mirrors top-level media files. |
| Property plots still show FastSpecFit SFRs | Repoint `WEDGE_PARQUET` to the CIGALE-HZ parquet from `build_cigale_rejoin.py`. |
| `plot_env_mass_continuous.py` cannot import `plot_style` | Script hard-codes the NERSC Illustris home path; other plot scripts use `ILLUSTRIS_ROOT`. |
| Gate G1 assert fails | SI self-eval npz test rows / truth must match the SI training cache masks. |
| Slurm graph products do not match `JRAPH_INPUTS_expanded_wedge.txt` | `sbatch_desi_bgs_bright_pipeline.sh` still targets the older narrow wedge dirs/sky cut; retarget or use the expanded-wedge runbook commands. |
| Cluster-deficit diagnostic scripts read the wrong run | Many FoG/coverage/shape helpers hard-code baseline `desi_wedge_flowjax_linear` paths — retarget before comparing to SI production. |
| `GalaxyCatalog` / `galaxy_catalog.py` KeyError on `RA` or `LOGMSTAR` | You passed the maglim zall FITS. That path uses `TARGET_RA`/`TARGET_DEC` and has no stellar mass. Use `loa-combined-lowz.fits`. |
| Gudhi `ModuleNotFoundError` on compute after login-node success | `srun` did not isolate PYTHONPATH / did not use `.../cosmic_env/bin/python`. |
| Wedge GNN arrays do not match `JRAPH_INPUTS_expanded_wedge.txt` | cuGraph default prefix is `desi_delaunay_cugraph`; validated products use `--out-prefix desi_bgs_cugraph`. |
| Maglim rebuild read a different zall | `build_bgs_maglim_catalog.py` defaults to the public DR2 zall; `load_catalog.py` uses `DESI_ZCAT_FILE` (collaboration redux). |
| `load_catalog.py` dies on a missing FastSpecFit file | Input list is hard-coded (`FASTSPEC_CATALOGS`), not a glob of `--fastspec-path`. |
| P12-A preflight exits 2 after hashes match | Expected while handoff/release gates are open. See `docs/p12a_vac_operations.md`. |
| P12-A worker rejects the git tree you just edited | `run_production.sh` executes `$ROOT/source/...` and checks `SOURCE_MANIFEST.json`. |
| `from production_inference import *` fails | Run `python workflows/p12a_vac/production.py`. `python -m` does not put that directory on `sys.path`. |
| Full-survey FITS used as a release catalogue | Header `PROVIS`/`SELSHIFT`; quality bit 32 on every row; bit 64 for `Z>=0.35`. `science_release_ready` is false. |
| P12-A class names disagree with GAT | Receipt `class_order` is void/sheet/filament/knot. GAT columns say wall/cluster for indices 1 and 3. |
| Shard JSON schema says `canary-vac-v1` | Production parts reuse that schema string. The full-survey marker is `FULL_VAC_COMPLETE.json`. |
| Loa source audit refuses on a login node | `p12a_loa_source_audit.py` requires `SLURM_JOB_ID` and a `nid` hostname. |

## Notes

- Active workflows and statuses: `ACTIVE_WORKFLOWS.md`
- SBI workflow depth: `workflows/sbi_inference/README.md`
- Reorg notes: `docs/migration/WORKFLOW_REORG.md`
- Shim policy: `scripts/SHIM_DEPRECATION.md`
- Prefer canonical `workflows/...` paths in new docs and launch scripts.
