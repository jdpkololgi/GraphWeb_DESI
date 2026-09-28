# P12-A Loa VAC — developer operations

Current application path for frozen Illustris U-PATCH + P12-A FMPE posteriors on
DESI BGS Loa. This page is the interface and pitfall index. Dated science
write-ups stay in the notes linked below; they are evidence, not launchers.

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
| Alignment screens | `workflows/p12a_vac/alignment_*.py` and matching `plot_alignment_*.py` |

`p12a_golden_replay.py` imports Illustris from the hard-coded home path
`/global/u2/d/dkololgi/TNG/Illustris`, not `ILLUSTRIS_ROOT`.

## Alignment investigation (do not rerun as a VAC)

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
