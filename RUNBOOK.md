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

Canonical:

```bash
python workflows/catalog/load_catalog.py --help
```

Compatibility wrapper:

```bash
python load_catalog.py --help
```

Example run (low-z GAT product, `z ≈ 0.01–0.06` only):

```bash
python workflows/catalog/load_catalog.py \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/catalogs
```

`--out-dir` defaults to `GRAPHWEB_CATALOG_DIR` (pscratch). Do not write these
FITS products under the repo or `$HOME` — home previously hit its 40 GiB quota.
Resolve paths via `shared/config_paths.py`:

- `GRAPHWEB_CATALOG_DIR`
- `GRAPHWEB_CATALOG_PATH`
- `GRAPHWEB_CATALOG_ZFLAGS_PATH`
- `GRAPHWEB_CATALOG_FASTSPEC_PATH`

This low-z catalog does **not** overlap the GraphWeb-BGS VAC redshift range
(`0.15–0.55`) and cannot supply `ABSMAG` for that work.

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
python workflows/graph_inference/graph_catalog.py --cache-mode cache-only --no-summary-plot
```

If model checkpoint path differs from config default:

```bash
python workflows/graph_inference/graph_catalog.py --model-path /path/to/trained_gat_model.pth
```

### Gudhi/cuGraph/Jraph wedge inference

This path is for DESI wedge eigenvalue regression with an Abacus-trained Jraph
model. It is not a replacement command for the all-sky GAT VAC workflow above.
Run it on Perlmutter from an environment with DESI, Gudhi/cuGraph where needed,
JAX/Haiku/Jraph, and access to the Illustris training repo.

Build the bright BGS catalog if it does not already exist:

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
`--scale-invariant-features`. For transfer diagnostics only, the script also
exposes `--edge-domain-adapt` and `--node-domain-adapt`.

Generate DESI diagnostics after inference:

```bash
RUN_DIR=/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_flowjax_linear_si

python workflows/sbi_inference/plot_desi_wedge_flowjax.py \
  --preds-npz "${RUN_DIR}/desi_wedge_flowjax_preds.npz" \
  --summary-json "${RUN_DIR}/summary.json" \
  --output-dir "${RUN_DIR}" \
  --no-regression
```

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

## Troubleshooting

| Symptom | Likely cause / fix |
| --- | --- |
| Home disk full / cannot commit | Catalogs and large outputs belong on pscratch (`GRAPHWEB_CATALOG_DIR`, `GRAPHWEB_SCRATCH_ROOT`). |
| FlowJAX abort on coordinate units | Rebuild DESI wedge with `--coord-units mpc`; do not use legacy Mpc/h products. |
| SI vs baseline mismatch | Match `--scale-invariant-features` to the training cache/model; the checked-in launcher is baseline-only. |
| Edge scaler looks wrong | Pass path1 fiberassign Abacus arrays to `--abacus-gnn-arrays`, not the regression wedge. |
| `desi_absmag_kcorr.py` `pkg_resources` ImportError | Use the repo script (it shims `pkg_resources`); or pin setuptools `<81` only if you must call DESI_ke outside this wrapper. |
| `add_ke` never finishes in interactive QOS | Keep the default `--max-rows 400000` z-cut subsample; full-catalogue k+e is a separate long job. |
| Figures missing after sync | `sync_figures_to_canonical.sh` only mirrors top-level media files. |
| Property plots still show FastSpecFit SFRs | Repoint `WEDGE_PARQUET` to the CIGALE-HZ parquet from `build_cigale_rejoin.py`. |
| `plot_env_mass_continuous.py` cannot import `plot_style` | Script hard-codes the NERSC Illustris home path; other plot scripts use `ILLUSTRIS_ROOT`. |
| Gate G1 assert fails | SI self-eval npz test rows / truth must match the SI training cache masks. |

## Notes

- Active workflows and statuses: `ACTIVE_WORKFLOWS.md`
- SBI workflow depth: `workflows/sbi_inference/README.md`
- Reorg notes: `docs/migration/WORKFLOW_REORG.md`
- Shim policy: `scripts/SHIM_DEPRECATION.md`
- Prefer canonical `workflows/...` paths in new docs and launch scripts.
