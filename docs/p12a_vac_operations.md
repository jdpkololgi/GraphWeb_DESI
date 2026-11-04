# P12-A Loa VAC — developer operations

Current application path for frozen Illustris U-PATCH + P12-A FMPE posteriors on
DESI BGS Loa. This page is the interface and pitfall index. Dated science
write-ups stay in the notes linked below; they are evidence, not launchers.
Mock repair is stopped. Later property and field figures do not replace the
frozen FITS.

Science plan (Illustris): `../../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md`.
Population-alignment subplan: `plan_mock_to_loa_alignment_investigation_20260926.md`.
Chronological handoff log: `p12a_vac_execution.md`.

## What exists, and what it is not

| Product | Code | Root on pscratch | Status in code |
|---|---|---|---|
| Frozen-artifact preflight | `workflows/catalog/p12a_vac_preflight.py` | report JSON only | Exit 2 while handoff gates are open. Does not open catalogues or checkpoints. |
| 16-core canary | `workflows/p12a_vac/loa_trial.py` | `.../outputs/p12a_loa_canary_20260925_v1` | `science_release_ready: false`. Geometry-selected cores, not a random sample. |
| Full-survey provisional VAC | `workflows/p12a_vac/production.py` | `.../outputs/p12a_loa_full_20260925_v1` | Technical merge marker `FULL_VAC_COMPLETE.json`. Header `PROVIS=True`, `SELSHIFT=True`. Not a science release. |

`docs/p12a_full_production_20260925.md` records 5,436,413 unique TARGETIDs at
`0.15 <= z < 0.55` (128 parts, 20-voxel cores, 48-voxel halo context,
alignment 8, 512 joint draws per supported galaxy). The planner itself stores
`ready['active_rows']` rather than that literal. Draws are per-galaxy joint
eigenvalue samples. They are not a joint posterior over the survey.

Do not treat GAT (`graph_catalog.py`), Jraph wedge, or FlowJAX wedge products
as this VAC. Those graphs are not the P12 count/response fields.

## Selection contract

`workflows/catalog/p12a_observation_patch.py` `successful_rows` (real data):

- finite `Z_not4clus > 0`
- `ZWARN == 0`
- `DELTACHI2 >= 25`
- `SPECTYPE == GALAXY` after strip

`loa_trial.prepare` then keeps context rows with `0.1 <= z < 0.6`, excluding
`0.585 <= z < 0.595`. Output shells are `0.15 <= z < 0.55` (`shell >= 0`).
Unit counts only: no PIP, ZFAIL, or SYS weights.

This is a third catalogue, distinct from:

- GAT low-z FastSpecFit (`load_catalog.py`): `0.01 <= z <= 0.06`, columns `RA`/`DEC`.
- Gudhi/Jraph/SBI maglim (`build_bgs_maglim_catalog.py`): no redshift cut,
  BGS_BRIGHT bits, columns `TARGET_RA`/`TARGET_DEC`.

The installed LSS success helper uses `DELTACHI2 > 40`. The `>= 25` rule is
the pinned GraphWeb cut (`p12a_quality_cut_audit.py` censuses both). Do not
swap them silently.

Coordinates use the frozen 17,001-point Planck18 lookup in `observer_xyz`
(`z <= 0.85`). Direct `comoving_distance` integration is not a drop-in
replacement. Cells are 5 Mpc.

## Canary CLI

Run the file path so the script directory is on `sys.path` (sibling imports).
`ILLUSTRIS_ROOT` overrides the default Illustris checkout. Candidate weights
are hard-coded under
`/pscratch/sd/d/dkololgi/abacus/p10_multiphase/p12a_halo48_candidate_20260924_v1`.

```bash
PY=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
ROOT=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1
$PY workflows/p12a_vac/loa_trial.py prepare  --root "$ROOT"
$PY workflows/p12a_vac/loa_trial.py finalize --root "$ROOT"
$PY workflows/p12a_vac/loa_trial.py infer    --root "$ROOT"
$PY workflows/p12a_vac/validate_canary.py "$ROOT"
```

`prepare` raises `FileExistsError` if `INPUTS_READY.json` already exists.
`infer` writes `DESI_LOA_P12A_HALO48_CANARY_VAC.fits` with `clobber=False`.
Canary core choice is truth-free: one median-occupancy core and one
smallest-boundary core per cap and shell, each with at least 16 galaxies.

## Full-survey CLI

`production.py` reads inputs only from the canary root above. `plan` refuses
if `PLAN.json` exists. Workers and merge call `check_sources` against
`SOURCE_MANIFEST.json` inside the **run root**, not the git checkout.

```bash
PY=/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python
ROOT=/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1
$PY workflows/p12a_vac/production.py plan      --root "$ROOT"
$PY workflows/p12a_vac/production.py benchmark --root "$ROOT"
$PY workflows/p12a_vac/production.py worker    --root "$ROOT" --part "$SLURM_ARRAY_TASK_ID"
$PY workflows/p12a_vac/production.py merge     --root "$ROOT"
```

Slurm wrapper `workflows/p12a_vac/run_production.sh` is not a login-node
shortcut. It requires a production root whose `source/` tree matches the
frozen manifest, sets `ILLUSTRIS_ROOT` to that frozen candidate snapshot,
unsets `PYTHONPATH`, and runs `$ROOT/source/workflows/p12a_vac/production.py`.
Editing this repository does not change an already frozen run. GPU: one A100
(`gpu&hbm80g`), `shared` QOS, 32 CPUs, 48 GB, 45 minutes. Stage `worker`
needs `SLURM_ARRAY_TASK_ID`; any other stage argument runs `merge`.

A completed core is reused only after input, ownership, and output hashes
match. A truncated completed file is a hard error. Seed is
`20260925 + global core_id`. `matmul` TF32 is off; cuDNN TF32 is on.

`check_benchmark.py` is the independent resume/parity check.
`plot_posteriors.py ROOT OUTDIR` is descriptive only: no DESI truth or
calibration claim.

## FITS columns and quality bits

Extension `ENV_POSTERIOR`. Header: `MODEL=P12A_H48`, `NSAMPLE=512`,
`LTHRESH=0.2`, `RSMOOTH=7`, `RUNIT=Mpc/h`, `ZTARGET=0.2`, `PROVIS=True`.
The merged full survey also sets `SELSHIFT=True`.

Supported rows store mean and quantiles 5/16/50/84/95 of the 512 ordered
eigenvalue draws, plus `P_WEB` (4,). `P_WEB[:, j]` is the fraction of draws
whose count of eigenvalues above 0.2 equals `j`. Canary receipt class order
is `void, sheet, filament, knot`. The older GAT pickle uses the words
wall and cluster for indices 1 and 3. Same index convention, different names.
Unsupported rows keep `TARGETID` and null summaries.

`QUALITY` bits written at inference (`loa_trial.py` / `production_inference.py`):

| Bit | Meaning |
|---:|---|
| 1 | Unsupported; posterior null |
| 2 | Sparse shell `z >= 0.45` (`shell == 3`) |
| 4 | Boundary distance `< 7/0.6766` Mpc (one smoothing length) |
| 8 | Boundary distance `< 2 * 7/0.6766` Mpc |
| 16 | A posterior feature outside the training min/max envelope |
| 32 | Provisional. Set on every row. |

Merge then ORs bit 64 when `Z >= 0.35`: unresolved selection-density excess,
not a correction and not marginalized. Shard FITS written before merge do
not have bit 64. Completion JSON is the bit legend for that product.

## Script map

Invoke as file paths. Most alignment and closure scripts take `--root` and
write evidence JSON; they do not update the VAC.

| Role | Entry |
|---|---|
| Small-artifact preflight (no FITS) | `workflows/catalog/p12a_vac_preflight.py` |
| Observer lattice / quality cut | `workflows/catalog/p12a_observation_patch.py` |
| Loa hash recheck (needs a `nid` Slurm job) | `workflows/catalog/p12a_loa_source_audit.py` |
| Full data + 18 randoms hash (needs `SLURM_JOB_ID`) | `workflows/catalog/p12a_refreeze_full_sources.py` |
| Quality-cut census, no selection change | `workflows/catalog/p12a_quality_cut_audit.py` |
| Golden mock replay, no DESI sampling | `workflows/catalog/p12a_golden_replay.py` |
| Canary / full survey | `loa_trial.py`, `production.py` |
| Read-only closure | `workflows/p12a_vac/closure_*.py`, `closure_summary.py` |
| Exposed-ph006 selection experiment | `mock_selection_test.py` (`prepare`/`sweep`), `run_mock_test_phase.sh` |
| Alignment screens (stopped) | `workflows/p12a_vac/alignment_*.py` and matching `plot_alignment_*.py` |
| Stopped v1-forward / selection helpers | See [below](#stopped-forward-and-selection-helpers). Do not invoke `prepare_v1_forward.py` with its defaults. |
| Spatial atlas of the frozen FITS | `workflows/p12a_vac/visual_atlas.py` |
| Wedge CIGALE cache vs current VAC | `workflows/p12a_vac/property_environment.py` |
| Full-survey CIGALE join | `loa_cigale_full.py`, `loa_cigale_product.py`, `loa_cigale_plots.py` |
| Matched P12-A vs CFM properties | `workflows/p12a_vac/field_property_environment.py` |
| Offline Loa density viewer | `loa_field_viewer.py`, `check_loa_viewer.py` |

`p12a_golden_replay.py` imports Illustris from the hard-coded home path
`/global/u2/d/dkololgi/TNG/Illustris`, not `ILLUSTRIS_ROOT`.

## Alignment investigation (stopped; do not rerun as a VAC)

Stopped by user direction on 2026-09-28. The plan is retained evidence, not an
execution queue. Do not launch further mock repair, population switches, or
retraining for this provisional VAC.

Canonical plan and latest checkpoint:
`plan_mock_to_loa_alignment_investigation_20260926.md`.
Independent controls (boundary erosion, DR1 vs Loa, velocity convention):
`p12a_independent_controls_20260926.md`.
Producer N/S and Holi follow-up: `p12a_producer_followup_20260926.md`.

Constraints encoded in those scripts:

- Read-only. No mock regeneration, new phase, retraining, or VAC replacement.
- Producer N and S files overlap in halo location and are not galactic caps.
  Do not concatenate them.
- Numerical `r < 19.5` cuts are not proof of passband or `BGS_TARGET` parity.
  The inspected Y3 forFA `BGS_TARGET` columns are zero.
- Holi close n(z) can be a construction target (thinning to Loa n(z)), not
  independent luminosity-model agreement.
- Public dipole bright/faint splits are not DESI `BGS_BRIGHT` / `BGS_FAINT`.

## Stopped forward and selection helpers

These files are the 2026-09-25–28 mock-repair diagnostics. Preparation job
59013807 was cancelled with that queue. None of them rewrite the frozen halo48
FITS. The measurements stay in the dated notes; this section is how the
scripts behave if they are opened.

Read, do not relaunch: `p12a_loa_completeness_forward_20260928.md`,
`p12a_v1_population_diagnosis_20260928.md`,
`p12a_joint_colour_luminosity_20260928.md`,
`p12a_canonical_v1_provenance_20260928.md`,
`p12a_published_lf_20260928.md`, `p12a_dr2_recipe_review_20260928.md`,
`p12a_selection_diagnostics_20260925.md`,
`p12a_selection_test_protocol_20260925.md`.

### `prepare_v1_forward.py` defaults to the DA2 mock tree

This is a vendored LSS bright-mock adapter. The `NERSC_HOST` check and
`ArgumentParser` sit at module level, so importing the file parses `sys.argv`.
Off Perlmutter or Cori it calls `sys.exit`. `NERSC_HOST` must be `cori` or
`perlmutter`.

Defaults that matter:

- `--isProduction` defaults to `y`. That branch forces `base_output` to
  `/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks` and sets overwrite false.
- `--base_output` itself defaults to that same CFS directory.
- `--downsampling` defaults to `n`, `--apply_mask` to `n`, and `--prog` to
  `bright`. The realization loop is `range(realmin, realmax)`, so the
  defaults `0` and `1` build only realization 0. Both launchers pass those
  same bounds.
- `--overwrite` uses `type=bool`. Any non-empty string, including `False`,
  becomes `True`. The production branch overwrites that flag back to false
  after parsing. A non-production run does not. The table is written with
  `targets.write(..., overwrite=args.overwrite)`.

`process()` seeds NumPy with `20260928 + real`. Bright globals use
`specver='loa-v1'`.

The only in-repo launchers override the CFS default. Both
`prepare_v1_forward.sbatch` and `smoke_v1_forward.sh` pass
`--mockver ab_secondgen_cosmosim`, `--isProduction n`, `--apply_mask y`,
`--downsampling n`, `--rbandcut 19.5`, and `--base_output` equal to
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/v1_forward_ph000_20260928`.
The sbatch also runs `sha256sum -c` on
`docs/evidence/p12a_joint_colour_luminosity_20260928/FORWARD_SOURCE.sha256`
(those four hashes still match this tree) and asserts installed LSS `HEAD`
is `d942b990860e016e7558c740966fea3e5ce7e7f3`. It stops after preparation
plus `check_v1_preparation.py`. `smoke_v1_forward.sh` is the foreground
`salloc` command and continues through one altMTL action.

### Chain after a scratch forFA file

Smoke-script root: `$OUT/SecondGenMocks/AbacusSummitBGS_v2`, with
`$OUT` the ph000 scratch directory above.

1. `check_v1_preparation.py <forFA0.fits>` has no `__main__` guard. It
   writes `PREPARATION_READY.json` beside the FITS with
   `qualification: false` and `placeholder_zwarn: true`. It requires
   sequential `TARGETID` from 1, finite `SUBPRIORITY` in `[0, 1)`,
   `BGS_TARGET` exactly 2 (bright), 1 (faint), or 9 (faint HIP), bright
   `R_MAG_APP` below 19.54 where `DEC > 32.375` and below 19.5 otherwise,
   other rows in `[cut, 20.175]`, and `ZWARN == 0`. That ZWARN check is the
   upstream ledger placeholder, not a spectral-quality model.
2. `initialize_v1_forward.py <forFA> <outdir> BRIGHT` is positional and also
   runs at import. It refuses if `<outdir>/Univ000` already exists.
   `recreate_tileTrack` is false: it copies the DA2 BRIGHT tile tracker and
   refuses if any `DONEFLAG` is already set.
3. `run_v1_altmtl.py --root <that BGS_v2 directory> [--smoke]` asserts the
   root starts with the ph000 scratch prefix, requires
   `PREPARATION_READY.json` and the same LSS commit, and runs tile shells
   under `bash -e`. `--smoke` advances one tracker action and writes
   `SMOKE.json`. Without `--smoke` it writes `ALTMTL_COMPLETE.json` only
   when no actions remain.

### Other stopped writers

| Script | Writes | Constraint |
|---|---|---|
| `selection_audit.py` | `SELECTION_AUDIT_V2.json` under `p12a_loa_canary_20260925_v1` | No CLI. Counts cuts on the canary Loa file, DA2 mock6 BRIGHT, and the ph006 observed FITS. Ablations are not new science cuts. |
| `mock_assignment_audit.py` | `MOCK_ASSIGNMENT_AUDIT.json` in that canary root | **Runs on import.** Bins mock rows by parent redshift through `TARGETID`, including unassigned sentinels. |
| `selection_diagnostics.py` | `docs/figures/p12a_selection_20260925/` | Positional `counts` or `phase`. |
| `selection_sensitivity.py` | `.../outputs/p12a_selection_sensitivity_20260925_v1` (`exist_ok=False`) | GPU replay of eight exposed ph006 cores, four density arms, 512 draws. Imports `workflows.abacus_tweb` and `workflows.sbi`, which this repository does not contain. Not a DESI correction. |
| `input_distribution.py` | `<canary-root>/INPUT_DISTRIBUTIONS.json` | Positional canary root. The canary is geometry-selected. |
| `upstream_population_diagnostic.py` | `UPSTREAM.json` beside the selection figures | ph006 raw cutsky on the same common pixels. Magnitude loops are sensitivity tests. |
| `summarize_passband_followup.py` | `SED_ROBUSTNESS.json` and `VALIDATION.json` in `docs/evidence/p12a_passband_followup_20260926/` | Rewrites committed evidence and checks script hashes. |
| `verify_response_sampling.py --output` | The requested JSON | 512-row parity of Loa `response_at` against ph006 voxels. Raises if distance or support differs. Imports `workflows.sbi.p12_prepare_base_response_dataset`, absent here. |
| `verify_mock_test_baseline.py --root` | `BASELINE_REPLAY.json` | BRIGHT columns versus official `forFA6_nomask.fits`. FAINT RNG is not claimed. Seed recorded as 62026. |
| `join_mock_test_baseline.py --root` | `BASELINE_JOIN.json` | Exact RA/DEC/RSDZ join. Requires `SLURM_JOB_ID`. |
| `replay_mock_y5_footprint.py --root` | `Y5_REPLAY.json` | Y5 cut plus post-cut TARGETID reindex. Requires `SLURM_JOB_ID`. |
| `materialize_mock_y5_variant.py --root` | `y5_baseline/forFA6_nomask.fits` | Needs a passing `Y5_REPLAY.json`. `FileExistsError` if the FITS exists. Header `TESTONLY`. |
| `validate_mock_y5_product.py --root` | `MATERIALIZED_VALIDATION.json` | Read-back of 10,000 bright rows against the official forFA file. |

`plot_selection_sensitivity.py` and `plot_upstream_population.py` only draw
the JSON above into `docs/figures/p12a_selection_20260925/`.
`plot_mock_selection_test.py --root` and `plot_reference_mock_check.py --root`
write `selection_screen.png`, `FINE_BIN_REVIEW.json`,
`reference_comparison.png`, and `SUMMARY.json` beside the run directory they
are given.

## Notes after the frozen VAC

These notes do not change selection, weights, or posteriors. Read the dated
page for the measurement; this list is only the index.

| Note | Role |
|---|---|
| `p12a_provisional_release_20260928.md` | Release limits, FITS integrity, Planck18 **Mpc** spatial atlas (`visual_atlas.py`) |
| `p12a_mock_choice_20260928.md` | Provisional release kept; Abacus v1 is not a validated replacement |
| `p12a_historical_footprint_20260928.md` | Historical BRIGHT footprint does not explain the high-z deficit |
| `p12a_realised_lf_20260928.md` | Parent luminosity vs recovered LF |
| `p12a_cigale_environment_20260928.md` | Wedge CIGALE cache vs current `P_WEB` |
| `p12a_property_literature_assessment_20260928.md` | How to read those modest contrasts |
| `p12a_contrast_checks_20260928.md` | Frozen-model contrast recovery on saved ph006 and the wedge |
| `p12a_historical_model_comparison_20260928.md` | Pointer to the Illustris knot-model replay |
| `p12a_loa_cigale_20260929.md` | Full-survey CIGALE join and controlled contrasts |
| `p12a_transition_audit_20260929.md` | Pointer to the Illustris boundary-uncertainty audit |

Three CIGALE paths are different inputs:

- `property_environment.py` reads the pinned full-survey FITS
  (SHA256 `cf472cb4…f50ebfbc`) and the wedge parquet
  `desi_wedge_cigale_hz/desi_wedge_env_props.parquet`. `--test` is a synthetic
  check. A real run refuses another VAC hash and an existing `--out`.
- `loa_cigale_full.py` joins
  `/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits`
  by TARGETID, then drops pairs with separation ≥1 arcsec or
  `|dz|/(1+z)≥0.001`, and leaves repeated consistent rows unmatched
  (`CIGALE_AMBIGUOUS`). `loa_cigale_product.py` writes the enriched FITS.
  `loa_cigale_plots.py` plots that product. The original VAC file is hashed
  before and after.
- `field_property_environment.py` is the matched CFM comparison below. It does
  not read either of those plot caches.

## Matched CFM properties and the offline 3D viewer

`field_property_environment.py` compares frozen P12-A `P_WEB` with CFM field
marginals on identical galaxies. `loa_field_viewer.py` builds a self-contained
HTML view of the same atlas. Neither writes the VAC. Atlas tiles are
independent local conditionals. A shared sample index across tiles is not one
survey-wide posterior draw. Sky-block intervals are selected-sample
diagnostics. They omit SED covariance and model-transfer error.

Both catalogue entrypoints raise `RuntimeError` unless `SLURM_JOB_ID` is set.
`--out` is created with `exist_ok=False`. Provenance is enforced up front:
`ATLAS_COMPLETE.json` must have `passed`, `vac_properties.npz` must match
`PLAN.json` `properties_sha256`, and the field-marginal file must match the
atlas receipt. `TARGETID` values must be unique.

The comparison sample keeps a row only when it is supported, uniquely
CIGALE-matched, `CIGALE_SPECTYPE=GALAXY`, `field_eligible`, and both CG15 and
CG5 have finite positive mass and SFR, with CG15 log mass in `[8, 13)`.
Environment order is void, sheet, filament, knot. The low-redshift scope is
`0.20 ≤ z < 0.30`. Controlled means use common cells (redshift bins of 0.025
from 0.15, and 0.1 dex mass bins when mass is controlled) and one reference
mix. A cell is kept only when every model and class has effective
`N ≥ 20`. The 16–84% interval is 64 resamples of nside-8 sky blocks
(`seed 20260930`). `--details` is separate: raw mean log sSFR inside each cap's
owned core, 16 local draws, `0.20 ≤ z < 0.30`, pooled support `≥ 0.5`, and a
class interval only when every draw has at least 20 members. If the receipt's
minimum `core_support` is below 0.25, that cap is dropped. The λth panel
compares CFM `p_web` with `p_web_threshold0` on the same cells. It is not a
second P12-A model.

The viewer embeds Plotly and the arrays in `Loa_3D.html`. It does not fetch
remote scripts. Coordinates are the atlas equatorial Cartesian grid in
**Mpc/h**. That is not the Planck18 **Mpc** slab in `visual_atlas.py`. Display
cuts never change the analysis sample:

- Overview: 100,000 galaxies per cap, chosen by a stable TARGETID hash.
- Detail: 2× pooled cells of 13.532 Mpc/h, marching cubes at 1, 2, and 4 times
  the mean where pooled support is at least 0.5, and 30,000 galaxies.
- Atlas (omit `--details-only`): the same pooling, then a 1-in-4 sublattice of
  supported owned-core voxels. No smoothing across tile seams.

Density isosurfaces are a display threshold. They are not the λth = 0.2
T-web class used in the property plots. `--surface-library` is prepended to
`sys.path` for scikit-image. `check_loa_viewer.py` also requires
`SLURM_JOB_ID`. It fails the job on a page error, any `http(s)` request, or
horizontal overflow at 430×932.

Perlmutter launcher (CPU, 32 cores, 30 minutes). The first argument is an
immutable source snapshot that contains `SOURCE_SHA256.json`. The script checks
those hashes, then requires atlas `rows=5436413`, `tiles=2472`, and `draws=8`.
Editing the git checkout does not change that snapshot.

```bash
ATLAS=/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/loa_field_atlas_20260930_v1
sbatch workflows/p12a_vac/submit_field_properties_20260930.sbatch /path/to/immutable/snapshot
```

That sbatch writes `$ATLAS/property_comparison_v1` and `$ATLAS/viewer_full_v1`
(atlas views included). It does not pass `--details`. Local draw variation is
a second command:

```bash
python workflows/p12a_vac/field_property_environment.py \
  --root "$ATLAS" --out "$ATLAS/local_draws" --details
```

`--test` on `field_property_environment.py` checks the paired estimator only.
`test_field_property_environment.py` and `test_loa_field_viewer.py` need a
compute allocation. The viewer smoke uses three real smoke tiles, writes a
disposable HTML, and records `SMOKE_SCOPE.json` so that fixture is not a
science atlas.

Committed evidence in `docs/figures/loa_field_properties_20260930/` is the
**details-only** viewer plus the local-draw figure. `VIEWER.json` has NGC and
SGC detail and `global_joint_posterior: false`. It has no atlas views. The
browser check passed with no remote requests and no mobile overflow. NGC shows
30,000 of 310,030 galaxies (18,225/18,432 supported display cells); SGC shows
30,000 of 170,611 (12,198/18,432). `legacy_density_replay` is true for NGC and
false for SGC. The HTML (34,255,787 bytes in that receipt) is not in git.
`LOCAL_RESULTS.json` records raw, not mass/redshift-controlled, log sSFR for
12,074 NGC and 9,432 SGC owned-core galaxies. The matched eight-figure
`RESULTS.json` from `main()` is not in this checkout.

## Pitfalls

| Symptom | Cause |
|---|---|
| Preflight exit 2 after hashes match | Expected while release gates are open. Other exceptions are execution failures. Reports are never overwritten. |
| `FileExistsError` on prepare/plan/audit output | Those writers refuse to replace an existing marker or report. |
| Worker hash failure after a git edit | The Slurm script runs the frozen `$ROOT/source` tree. |
| `import production_inference` fails | Run `python workflows/p12a_vac/production.py`, not `python -m`. |
| Canary FITS not replaced | `fitsio.write(..., clobber=False)`. |
| Bit 64 missing on a shard | Only `production.merge` sets it, on `Z >= 0.35`. |
| Class names disagree with GAT plots | P12-A receipt says sheet/knot; GAT columns say wall/cluster. |
| Full survey used as a paper catalogue | `science_release_ready` is false. Bit 32 is universal. Selection shift is flagged, not corrected. |
| Login-node catalogue hash | `p12a_loa_source_audit.py` requires `SLURM_JOB_ID` and a `nid` host. |
| Mixing this VAC with wedge npz | Different selection, features, and model. |
| Part `VAC_COMPLETE.json` schema says canary | `production_inference.py` reuses schema `desi-loa-p12a-halo48-canary-vac-v1` for shards. Read `scope` and the merged `FULL_VAC_COMPLETE.json` (`schema` `p12a-loa-full-survey-v1`). Part files are `VAC_SHARD.fits` and `shards/core_{id:06d}.*`. The canary writer uses `core_{id:02d}_draws.npz` and `DESI_LOA_P12A_HALO48_CANARY_VAC.fits`. |
| `field_property_environment.py` or `loa_field_viewer.py` raises `RuntimeError` | Those entrypoints require `SLURM_JOB_ID`. `--test` on the property script does not. |
| `FileExistsError` on a property or viewer `--out` | The directory must not already exist. |
| Atlas sample index read as one cosmic draw | Tiles are independent. `VIEWER.json` sets `global_joint_posterior` false. |
| Viewer box overlaid on the release atlas | Viewer axes are Mpc/h. `visual_atlas.py` uses Planck18 Mpc. |
| Density surface read as a T-web class | Isosurfaces are a display cut. Classes in the property script use λth = 0.2. |
| CFM λth = 0 panel treated as P12-A | That panel compares two CFM definitions on the same cells. |
| sbatch ignored a git edit | It runs the snapshot passed as `$1` after checking `SOURCE_SHA256.json`. |
| Committed `VIEWER.json` described as the full atlas | The file in git has detail views only. The sbatch's `viewer_full_v1` is the atlas HTML, and it is not committed. |
| Bare `python workflows/p12a_vac/prepare_v1_forward.py` | Defaults are `--isProduction y` and the DA2 mocks CFS directory. The sbatch and `smoke_v1_forward.sh` pass `--isProduction n` and the ph000 scratch root. |
| `--overwrite False` on that preparer | `type=bool` turns any non-empty string into `True`. Production mode forces overwrite false afterwards; a scratch run does not. |
| Import of `prepare_v1_forward.py`, `initialize_v1_forward.py`, `check_v1_preparation.py`, or `mock_assignment_audit.py` | All four execute at import. The first two read `sys.argv`. The assignment audit writes canary JSON. |
| `run_v1_altmtl.py` raises on LSS `HEAD` or `--root` | It requires LSS `d942b990860e016e7558c740966fea3e5ce7e7f3`, `PREPARATION_READY.json`, and a root under `v1_forward_ph000_20260928`. |
| `selection_sensitivity.py` or `verify_response_sampling.py` cannot import | `workflows.abacus_tweb` and `workflows.sbi` are not modules in this repository. |
| `summarize_passband_followup.py` dirties git | It rewrites `VALIDATION.json` under `docs/evidence/p12a_passband_followup_20260926/`. |
