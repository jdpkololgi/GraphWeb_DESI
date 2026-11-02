# Six-fit Loa field operations

Exploratory observed-Loa density fields from six frozen fits. This is not the
P12-A per-galaxy VAC, not a retraining, and not a calibrated survey posterior.
Every production receipt in this chain sets `posterior_qualified` to false.

Two products share the same checkpoints and the prepared tile bundles, and they
must not be described as the same draw:

| Product | What is sampled | Noise namespace |
| --- | --- | --- |
| Independent-tile survey | One 32×16×16 core per atlas tile | `loa-survey-v1` plus tile id, seed, draw, stage |
| Shared-state cap | One synchronous latent canvas per cap and seed | `loa-shared-sixfit-v1` plus cap, seed, stage, draw index |

Tile identity is in the survey seed and absent from the shared-state seed.
Worker, shard, and requested draw count are in neither. The survey namespace
exists so a cap-only seed is not copied across 2,472 tiles. The shared-state
namespace exists so overlap consensus uses one noise field per cap.

Fits, in both products: legacy CFM (`legacy_v`), expanded U-Net endpoint
(`unet_x`), and transformer velocity (`transformer_v`), seeds 17 and 29.
Expanded checkpoints are step 9984 EMA. NFE 128 means 64 Heun steps. Coarse
draws are shared across the three families within a seed. The unit-density
exterior is the finite local-plus-wide operator, not a learned full-cap exterior.

The 2026-10-19 Abacus label-provenance warning in `CONTEXT.md` still applies.
Do not swap a corrected-label checkpoint into these bindings and call the
existing fields updated.

Campaign narratives, row counts, and Slurm accounting live in the dated notes
listed at the end. This page is the launch contract.

## Independent-tile survey

Root, hardcoded by `workflows/presentation/loa_survey_production_20261016.py`:

`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1`

Campaign id `loa-six-fit-20261016-v1`. The runner's own default `--root` is the
preflight directory `loa_survey_preflight_20261016_v2`, and its default
`--scope` is the 16-tile preflight. A production command must pass the
production root and an explicit tile list or shard. `initialize` refuses any
path other than the registered root and refuses an existing root.

Payloads require `SLURM_JOB_ID` and a hostname starting with `nid`
(`require_compute`). The runner does not itself count GPUs; the preview and
extension batch scripts request one 80 GB A100. None of these modules submit,
requeue, or retry a job.

```bash
# Lifecycle only. Does not sbatch.
python workflows/presentation/loa_survey_production_20261016.py initialize
python workflows/presentation/loa_survey_production_20261016.py validate-recovery
python workflows/presentation/loa_survey_production_20261016.py prepare
python workflows/presentation/loa_survey_production_20261016.py shard --shard 0 --draws 4

# Direct runner. --stop is exclusive. Default stop is 4; maximum is 32.
python workflows/presentation/loa_survey_runner_20261016.py infer \
  --root /pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1 \
  --tiles NGC_005_008_003 --first 0 --stop 4 --seconds 6600
```

`workflows/presentation/run_loa_survey_production_20261016.sh` is the batch
wrapper: it unsets `PYTHONPATH`, pins `cosmic_env`, and execs the production
module. `loa_survey_preview_20261016.slurm` is array `0-102%8`, draws 4, two
hours, `--no-requeue`. `loa_survey_extension_20261016.slurm` is the same array
at ten hours and draws 32. Its header says not to submit it until preview
review exists. The code agrees: `shard(..., draws=32)` raises `PermissionError`
unless `PREVIEW_REVIEW.json` has `extension_authorized` true and
`campaign_sha256` equal to the hash of `CAMPAIGN.json`. The preview audit
writes a different file, `PREVIEW_AUDIT.json`, and sets
`extension_authorized` false there. Passing the audit does not create
`PREVIEW_REVIEW.json`.

Index order is 103 shards of 24 tiles (2,472 tiles, including empty retained
sets). Draw bundles publish by renaming a hidden temporary directory onto a
path that must not already exist. Receipts use `os.link`, so an existing JSON
file is never replaced. A binding or checksum mismatch on an existing bundle
raises; it is not resampled. `--max-new` is an intentional stop for resume
tests. A wall-clock deadline keeps completed draws and does not retry.

Cleanup may unlink only non-`prepared/` members inside this campaign root, and
only after the archive's member hashes match. It does not recurse, and it does
not delete hidden partial directories. `pack_loa_survey_tile_20261016.py`
refuses an existing archive. Restore fills missing members through temporary
files and rejects links and path escapes. A second restore must preserve
modification times.

Parent selection used here is the half-cell interpolation edge cut on the
2,472-tile atlas (retained rows are not the full 5,436,413 VAC). Property
eligibility is a further mask. Report both; do not describe either product as
inference for every parent galaxy.

## Shared-state caps

Engineering and production entry points:

- `workflows/presentation/loa_shared_sixfit_20261026.py` — `smoke` or `benchmark` only.
- `workflows/presentation/loa_shared_field_io_20261026.py` — HDF5 writer and coverage checks. No fitting.
- `workflows/presentation/freeze_loa_shared_production_20261026.py` — hash-bound snapshot. Never submits.
- `workflows/presentation/loa_shared_production_20261026.py` — one approved cap and seed.
- `configs/loa_shared_production_20261026.json` — the approved task table.
- `workflows/presentation/loa_shared_production_20261026.slurm` — runs the unit tests, then production, from a frozen snapshot.

Smoke selects an 8-window geometric group (maximum summed row count; not chosen
from inferred values). Benchmark selects the full cap. Smoke uses NFE 128;
benchmark uses NFE 2. A benchmark file is a timing measurement, not a posterior.
Both modes require one GPU on an `nid` host, and `--out` must not exist.
JSON from this module is `open('x')`. Default `--seconds` is 3000; exceeding it
raises `TimeoutError` and leaves partial receipts.

```bash
python workflows/presentation/loa_shared_sixfit_20261026.py \
  --mode smoke --cap NGC --seed 17 --draws 1 --seconds 3000 \
  --out /pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/smoke_NEW
```

`freeze_loa_shared_production_20261026.py` refuses an existing
`source_production_v1` directory. It copies sources only when the GPU-tested
benchmark hashes still match the live files and the production sources match
`git show HEAD`. A dirty or uncommitted production file is a hard stop.
`ILLUSTRIS_ROOT` inside the batch script is the snapshot's `Illustris` tree,
not the mutable checkout.

Production reads that freeze, the approved JSON, and
`loa_survey_production_20261016_v1/INDEX.json`. It then creates

`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1/production_v1/<CAP>_seed<SEED>`

with `exist_ok=False`. Checkpoints and `PROGRESS.json` are written during the
run, and `PROGRESS.json` is the one receipt this writer replaces. They do not
enable resume: a second launch sees the directory and stops. About ten minutes
before the approved wall time the run raises `TimeoutError` and does not
resubmit. The approved table is four draws, NFE 128, one 80 GB GPU (the process
requires at least 70 GiB), TF32 off, expanded autocast bfloat16. Tasks are
NGC seeds at 32 hours and SGC seeds at 18 hours. The batch script accepts only
`NGC:17`, `NGC:29`, `SGC:17`, and `SGC:29`, and only the registered snapshot
path.

Expanded sampling must keep the native initial, proposal, and update
projections. A second velocity projection changes bfloat16 trajectories and
fails the single-window parity check (`max abs` must stay below `1e-6`).
`diagnose_loa_shared_parity_20261026.py` executes on import, needs one GPU and
`SLURM_JOB_ID`, and exclusively creates `PARITY_REPAIRED.json` (or
`PARITY_REPAIRED_UNCACHED.json`) under the shared-sixfit output root.

The HDF5 writer (`CELL = 6.766` Mpc/h) publishes with `os.link` and raises
`FileExistsError` if the field path exists. Uncovered cells are stored as NaN,
not filled with a mean. `posterior_qualified` is stored false. The coverage
preflight checks every parent galaxy in that cap from `vac_properties.npz`,
including rows the independent-tile product dropped at interpolation edges.
A shared-state coverage pass is not a statement that those dropped galaxies
have tile posteriors.

## Display packaging

Both packagers read the four-draw survey archives. They do not run inference
and they do not make the tiles jointly sampled.

- `pack_loa_fullcap_page_20261026.py` writes only
  `.../loa_survey_production_20261016_v1/page_fullcap_20261026_v2`, which must
  not exist. It requires `nid` and `SLURM_JOB_ID`. The display is the native
  central z plane (index 8), a four-draw mean, and 4×4 in-plane pooling, then
  16 log2 levels on [-4, 4]. Class counts are the saved draw counts, not
  labels decoded from the quantized plane.
- `export_loa_fullcap_site_20261026.py --output <new-dir>` is a different
  reduction: 2× pooling in all three axes, mean plus four draws, log2
  quantization from -8. The manifest states the mosaic is not a joint full-cap
  draw. Context galaxies are a deterministic retained-position sample, not
  environment-labelled targets.

`test_loa_fullcap_page_20261026.py` asserts `SLURM_JOB_ID` at import and serves
the existing page directory. It is not a login-node unit test.
`test_loa_survey_production_20261016.py` is a fixture test and does not read
the survey archives. `test_loa_shared_sixfit_20261026.py` and
`test_loa_shared_field_io_20261026.py` are synthetic; the batch script runs
them from the frozen snapshot before production.

## Where the dated notes stop

| Question | Read |
| --- | --- |
| Preflight budget, locks, 16-tile selection | `docs/loa_survey_preflight_20261016.md` |
| Preview array, archive recovery, extension gate | `docs/loa_survey_production_20261016.md` |
| Why four draws are not a joint posterior | `docs/e2e_fresh_eyes_briefing_20261016.md` |
| Two-tile property pilot | `docs/loa_model_page_pilot_20261015.md` |
| Approved shared-cap jobs and hashes | `docs/loa_shared_sixfit_20261026.json` and `configs/loa_shared_production_20261026.json` |
| Full-cap Page contents | `docs/loa_fullcap_page_20261026.md` |
| Tokyo slide scripts and the label warning | `docs/tokyo_source_inventory_20261027.md` |

## Filament-geometry follow-ups

These consume the frozen P12-A VAC. They do not change it. The early paragraphs
of `docs/filament_parent_stage_20261017.md` describe a plan; the execution
section at the end of that file is the completed parent extraction. Later
measurements are separate notes:

- Off-slab tracer shapes: `docs/filament_offslab_results_20261017.md`
- ph006 tidal truth versus strings: `docs/filament_mock_truth_20261017.md`
- Three matched examples: `docs/filament_mock_examples_20261017.md`
- Matter, sampling, and redshift space in those regions: `docs/filament_realspace_20261017.md`

`prepare_parent.py`, `compare_offslab.py`, `mock_truth.py`, `mock_examples.py`,
and `realspace_extract.py` require `SLURM_JOB_ID` and raise `FileExistsError`
if `--out` already exists. `realspace_compare.py` also requires `SLURM_JOB_ID`.
`realspace_sampling_null.py` runs on import and takes the results directory as
`sys.argv[1]`. `plot_offslab.py` overwrites figures under `--out`. ph006 is a
development mock that was already exposed; these scripts are not a blind phase.
