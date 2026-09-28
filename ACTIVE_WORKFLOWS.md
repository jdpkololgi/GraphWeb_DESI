# Active Workflow Index

## P12-A Loa application (updated 2026-09-28)

Current DESI application of the frozen Illustris U-PATCH + P12-A FMPE posterior.
Developer interface, FITS columns, quality bits, and pitfalls:
`docs/p12a_vac_operations.md`. Chronological handoff log:
`docs/p12a_vac_execution.md`. Science plan:
`../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md`.

- Preflight (no catalogue open): `workflows/catalog/p12a_vac_preflight.py`.
  Exit 2 means handoff gates are still open.
- Canary (16 geometry-selected cores): `workflows/p12a_vac/loa_trial.py`
  stages `prepare`, `finalize`, `infer`.
- Full-survey provisional VAC: `workflows/p12a_vac/production.py` stages
  `plan`, `benchmark`, `worker`, `merge`. Slurm wrapper
  `workflows/p12a_vac/run_production.sh` runs the frozen copy under the run
  root, not the live git tree.
- Both products set `science_release_ready: false`. Every row has quality
  bit 32. The merged catalogue also sets bit 64 for `Z >= 0.35` (unresolved
  mock–Loa count excess; not a correction).
- Active science follow-up is the mock-to-Loa alignment investigation
  (`docs/plan_mock_to_loa_alignment_investigation_20260926.md`), not a new
  inference entrypoint. Alignment scripts are read-only.

Legacy GAT/Gudhi/Jraph/FlowJAX workflows below remain usable for their
historical scope. They do not build P12 count/response fields.
Slurm, reserved phase-access and science-release gates remain in force.

This file is the Phase 0 quick reference for what to run in this repository now.
For scratch layout, conda envs, and Perlmutter job recipes, see `~/.claude/CLAUDE.md`;
for launch commands see `RUNBOOK.md`.

## Active

- Low-z catalog assembly for the GAT classifier:
  - Canonical: `workflows/catalog/load_catalog.py`
  - Compatibility shim: `load_catalog.py`
  - Writes to `GRAPHWEB_CATALOG_DIR` on pscratch by default (not the repo/home).
  - Hard-coded FastSpecFit filename list; `SPECTYPE==GALAXY` and `z ∈ [0.01, 0.06]`;
    columns are `RA`/`DEC`. Not a Gudhi/Jraph input.
- Bright BGS catalog assembly for Gudhi/Jraph/SBI:
  - Canonical: `workflows/catalog/build_bgs_maglim_catalog.py`
  - `ZWARN==0`, `DELTACHI2>=25`, `SPECTYPE==GALAXY`, BGS_BRIGHT bits; **no** z or
    mass cut. Columns are `TARGET_RA`/`TARGET_DEC`. Default `--zall-path` is the
    public DR2 zall, not `DESI_ZCAT_FILE`.
- GAT graph inference pipeline:
  - Canonical: `workflows/graph_inference/graph_catalog.py`
  - Compatibility shim: `graph_catalog.py`
  - Cache/VAC still default to the **repo**, not `GRAPHWEB_CANONICAL_CACHE_DIR`.
    Point `GRAPHWEB_CACHE_DIR` / `GRAPHWEB_VAC_OUTPUT_PATH` at pscratch before
    `rebuild`. `--no-summary-plot` still loads TNG300. See `RUNBOOK.md`.
- Gudhi/cuGraph graph construction:
  - Full graph: `workflows/graph_construction/build_desi_bgs_gudhi_graph.py`
  - Feature export: `workflows/graph_construction/desi_graph_features_cugraph.py`
  - Wedge subset: `workflows/graph_construction/subset_desi_graph_wedge.py`
  - Use Mpc coordinates (`--coord-units mpc`, default) for Abacus/Jraph/SBI parity.
  - Hemisphere-split AlphaComplex: **no north–south edges**. Keep maglim row
    order through graph → features → wedge. Validated cuGraph prefix is
    `desi_bgs_cugraph` (script default is `desi_delaunay_cugraph`).
  - Perlmutter Slurm wrapper: `workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh`
    (`STEP=2|2b|3`, or `SUBMIT_CHAIN=1`). Still hard-codes the older narrow
    wedge dirs/sky cut — see `RUNBOOK.md` before treating its outputs as the
    expanded Mpc-parity products.
- Jraph inference:
  - Active DESI wedge inference: `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`
  - Expanded-wedge launcher: `workflows/catalog/run_infer_expanded_wedge_mpc.sh`
  - Requires Abacus `node_feature_scaler` + 15-d increment reconstruction; see
    `RUNBOOK.md` / `workflows/catalog/BUILD_JRAPH_CACHES.md`.
  - Warns (does not abort) on Mpc/h metadata; defaults `JAX_PLATFORMS=cpu`;
    hparams `96/8/8/0.15` and `output_dim=15` are CLI/hardcoded, not in the
    checkpoint. Writes `desi_wedge_index_and_preds.npz`.
  - Legacy GAT-Delaunay pickle builder (not production):
    `workflows/jraph_inference/build_desi_wedge_jraph_cache.py`
- FlowJAX/SBI posterior inference:
  - DESI wedge NPE inference: `workflows/sbi_inference/infer_desi_wedge_flowjax.py`
  - Abacus self-inference reference: `workflows/sbi_inference/infer_abacus_self_flowjax.py`
  - DESI posterior diagnostics: `workflows/sbi_inference/plot_desi_wedge_flowjax.py`
  - Property-closure join/plots:
    `workflows/sbi_inference/build_desi_wedge_property_join.py`,
    `workflows/sbi_inference/plot_property_environment_closure.py`
  - Workflow README: `workflows/sbi_inference/README.md`
  - Requires the same Mpc-parity DESI wedge products, path1 Abacus edge-scaler
    arrays, and `ILLUSTRIS_ROOT` for the Illustris FlowJAX model code.
  - Aborts on Mpc/h metadata; uses GPU unless you force CPU; baseline mode
    asserts path1 edge-scaler constants. Default post-hoc λ sort; SI node
    medians skip Clustering (col 1).
- Utility:
  - Canonical: `workflows/utilities/galaxy_catalog.py` — GAT low-z helper
    (`RA`/`DEC`, `LOGMSTAR>=9`). `--help` is a stub. Do not point it at the maglim FITS.
  - Canonical: `workflows/utilities/investigate_edges.py` — GAT
    `DESI_alpha_geom.pt` edge-length QA, not Gudhi arrays.
  - Compatibility shims: `galaxy_catalog.py`, `investigate_edges.py`
- Shared helpers:
  - Config paths: `shared/config_paths.py` (shim: `config_paths.py`)
  - Abacus/Jraph coordinate and edge-feature parity: `shared/abacus_gnn_parity.py`
    (`mpc` vs legacy `mpc_per_h`; reverse edges negate direction and take
    `1/density_contrast`, then log+z-score cols 0 and 4). `--scale-invariant`
    divides `edge_length` by the per-graph median *before* the log.
    Inference scripts load this file by path so Illustris `shared/` does not
    shadow it.
  - Plot style: `shared/plot_style.py` / `PLOT_STYLE_GUIDE.md` (`GRAPHWEB_FONT_DIR`)
- NERSC helpers: `scripts/where_am_i.sh`, `scripts/check_quota.sh`,
  `scripts/desi_prods.sh`, `scripts/sync_figures_to_canonical.sh`,
  `templates/job.sbatch`. Slurm graph jobs need the env’s absolute python
  (`.cursor/rules/conda-env-srun-python-path.mdc`).

## Visualization

- Active Jraph path1 3D CWEB / embedding notebook:
  `workflows/visualization/visualize_desi_wedge_cweb_3d.ipynb`
  Edit cell 1 for paths. It reads `INFER_DIR/desi_wedge_index_and_preds.npz`
  (Jraph), not FlowJAX `desi_wedge_flowjax_preds.npz`. HTML is written under
  `INFER_DIR`. This is **not** the FlowJAX plotter
  (`plot_desi_wedge_flowjax.py`).
- FlowJAX talk visuals (Keynote): `build_skewer_animation.py` (real NPE pencil
  beam), `build_skewer_idealised.py` (synthetic, no catalogue),
  `render_skewer_video.py` (HTML → mp4/gif), `render_class3d_video.py`
  (FlowJAX `hard_class` fly-through GIF), `plot_eig_dist_vertical.py` /
  `plot_eig_dist_buildup.py`. Inverse of the notebook: these need FlowJAX
  preds, not the Jraph npz. See `workflows/sbi_inference/README.md`.
- `workflows/visualization/path1_desi_wedge_inference_summary.md` is a
  NERSC-only symlink into pscratch (`.../path1_epoch8056_on_desi_expanded_wedge_boxcox_fix/INFERENCE_SUMMARY.md`);
  it is dangling off Perlmutter and should not be treated as in-repo docs.
- Archived notebooks: `workflows/visualization/archive/` (LOA catalog, graph subvolume; see README)

## Experimental

- Notebook-driven variants under this repo are exploratory unless explicitly promoted.
- Candidate luminosity-channel support diagnostic:
  - `workflows/sbi_inference/desi_absmag_kcorr.py`
  - Computes DESI `ABSMAG_RP1` with the DESI LSS k+e implementation and prints
    a marginal comparison with a raw Abacus cut-sky sample.
  - This is not a production pass/fail gate: the current DESI and Abacus
    selections, redshift distributions, and magnitude conventions are not
    matched. See `RUNBOOK.md` before interpreting its output.
- Transfer / capacity gates (Abacus-domain; see `workflows/sbi_inference/README.md`):
  - `gate_g1_gnn_vs_gbm.py` — GNN vs HistGBM on identical SI features.
  - `gate_g15_g2_rsd_luminosity.py` — RSD penalty and luminosity-weight gain.
  - `measure_nz_mock_vs_desi.py` — path1 mock vs DESI BGS n(z) in a shared sky box.
- Property-science follow-ups (env vars `WEDGE_PARQUET` / `WEDGE_FIGDIR`):
  - `build_cigale_rejoin.py` — swap CIGALE-HZ mass/SFR onto SI wedge posteriors.
  - `plot_sfms_environment.py`, `plot_mstar_color_environment.py`,
    `plot_env_mass_continuous.py`, `investigate_env_property_signal.py`.
- Cluster-deficit falsification / transfer diagnostics (indexed in
  `workflows/sbi_inference/README.md`, plan in `CLUSTER_DEFICIT_FIX_PLAN.md`):
  FoG LOS alignment, training-support / SI coverage, shape misclassification,
  embedding MMD, λ_th sweep, edge-scale / recovery bars, skewer animations, plus
  velocity-dispersion / smoothing-scale / mass-anchored investigations. These
  are research diagnostics, not canonical launch paths; many hard-code baseline
  run directories.
- Historical Jraph unit-parity one-off (Mpc/h DESI arrays, 2026-05-29; do not
  re-run on current Mpc products):
  `workflows/jraph_inference/experiment_desi_feature_parity.py`

## Legacy/To retire

- Any duplicate notebook/script flow superseded by `graph_catalog.py`.
- Older Mpc/h DESI wedge artifacts are deprecated; prefer the expanded Mpc-parity
  paths in `workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt`.

## Compatibility policy

- Root-level script names are temporary wrappers to keep existing commands working.
- Prefer calling canonical workflow paths for new scripts, docs, and SLURM launchers.
- Config paths canonical module: `shared/config_paths.py` (shim: `config_paths.py`)
