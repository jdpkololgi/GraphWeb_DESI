# Loa six-fit full-redshift runner preflight

## Authorization and scope

The user approved one CPU-node-hour on `desi` and two GPU-hours on `desi_g`,
including release of both allocations, to implement/validate the runner and
report a measured budget **before any production submission**. No production
submission, training, checkpoint selection or Page publication is authorized
by this preflight. The existing Page and pilot products are unchanged.

The frozen fits are the earlier CFM, expanded U-Net clean endpoint and
transformer velocity, each seeds 17/29. Samplers retain NFE128 (64 Heun steps
for the earlier CFM), the common physical 32x16x16 cores and the finite
free-space local-plus-wide tidal operator, with unit density beyond the wide
domain. Independent tiles are not a coherent survey-wide posterior. Numerical
success is not posterior calibration.

## Implementation

- `workflows/presentation/loa_survey_runner_20261016.py`: one-time catalogue
  index/join, explicit tile/shard preparation, streamed six-fit inference,
  immutable source bindings and atomic per-draw output bundles.
- `workflows/presentation/run_loa_survey_preflight_20261016.sh`: sanitized
  environment/payload launcher; it contains no scheduler submission logic.
- `workflows/presentation/test_loa_survey_runner_20261016.py`: focused tests
  for seed addresses, shell edges, selection, all 103 shard boundaries,
  atomic commits, checksum rejection and exclusive/shared locks.
- `workflows/presentation/loa_survey_preflight_checks_20261016.py`: CPU
  concurrency/array-parity benchmark and fixed GPU stop/resume/extension/replay
  validation, followed by a cap/shell-weighted budget calculation.
- `workflows/presentation/loa_survey_empty_tile_check_20261016.py`: bounded
  one-tile check of inference with no retained galaxy rows.
- `workflows/presentation/pack_loa_survey_tile_20261016.py`: atomic per-tile
  archive export with verification of every member checksum, retaining all
  input bundles. These archives contain private galaxy data, not Page payloads.

The root manifest binds the runner, scientific source tree, frozen checkpoints,
atlas plan, observed response inputs and catalogue. Each prepared tile is
checksummed. Each draw binds its tile, seed, model, draw index and shared coarse
array. A directory rename commits arrays and receipt together. Completed
bundles are validated before reuse; mismatches fail closed. Hidden incomplete
temporary directories are ignored, not deleted or treated as complete.

Random addresses include the **tile ID**, seed, draw and stage. They exclude
worker, shard and requested draw count. This deliberately replaces the old
two-region pilot's cap-only addressing, which would repeat noise across tiles
if copied unchanged to the survey. Expanded families share fine noise, and all
families share coarse draws within each seed; legacy fine noise has a different
shape and its own address. This does not make tiles a joint realization.

No script submits/requeues jobs or automatically retries failed work. Existing
valid draws are reusable only with identical bindings. Future source edits must
not be silently accepted against these receipts.

## Coverage and selection

The parent manifest has 2,472 tiles: 1,626 NGC and 846 SGC, across supported
redshift 0.15-0.55. The full unique-TARGETID photometry join checks exact RA,
Dec and redshift agreement and requires one match per parent row. Duplicates
within a single input chunk are counted correctly with indexed accumulation.

Of 5,436,413 parent rows, 4,628,192 survive the existing half-cell interpolation
edge rule; 808,221 are excluded before support/property cuts. There are 49 tiles
with no retained galaxy rows. Do not describe this as inference for every parent
galaxy or fill those gaps by unregistered interpolation/blending.

The 16-tile preflight selects two tiles per cap and median-redshift shell
([.15,.25), [.25,.35), [.35,.45), [.45,.55]). Within each cap/shell, edge/interior
are the lower/upper quartiles of median catalogue BOUNDARY_MPC. This is a
relative geometry stratification, not a claim of pristine interior support.
Selection uses no inferred fields or CIGALE outcomes. Some selected high-z
tiles have very few galaxies and one has no property-eligible galaxies; these
cases are retained and reported rather than outcome-filtered.

Property eligibility retains unique unambiguous CIGALE GALAXY matches, saved
catalogue support and interpolated support >=0.5, 9<=log10(Mstar/Msun)<12 and
corrected r<19.5. North photometry uses the registered -0.04 mag correction.
Property-specific validity remains necessary during aggregation; no missing
property value is implicitly counted as quenched.

## CPU validation and retained failure

CPU allocation `59833970`, node `nid200006`, used `desi`, `interactive`, one
node, requested 64 logical CPUs and one hour, with Scratch access. It was
released after 11m21s; allocation state `COMPLETED`, exit 0:0.

The initial serial preparation passed, but its concurrency test exposed an
exclusive catalogue-index lock being used by read-only workers. Seven workers
correctly failed instead of racing. The creation lock remains exclusive;
completed-index readers now take a shared lock, and tile writers remain
exclusive. A regression test covers this distinction. All nine focused tests
pass after the fix.

The failed v1 attempt is retained, not overwritten. Fresh v2 preparation binds
the corrected source. Its index took 188.883 s, including a 171.307 s full
photometry scan/join. All 16 expanded observations reproduce the legacy central
observations exactly (maximum absolute difference zero).

Serial preparation of 16 tiles took 46.377 s. Eight four-thread workers,
preparing the same 16 tiles into a separate benchmark directory, took 14.573 s
including process startup and index checks: measured speedup 3.182x. Every saved
array is exactly equal to the serial result. This is a small concurrency
benchmark, not a sustained full-survey filesystem guarantee.

## GPU validation and measured budget

GPU allocation `59834307`, node `nid200273`, used `desi_g`,
`shared_interactive`, one NVIDIA A100-SXM4-80GB, 32 logical CPUs, a two-hour
limit and Scratch access. Both GPU payload steps completed with exit 0:0.
It was released after **1h17m37s**, allocation state `COMPLETED`, exit 0:0.
Final queue inspection was empty. Actual allocation use was 0.18917 CPU-node-h
and 1.29361 GPU-h, within the approved one/two-hour limits. Production has not
been submitted.

The fixed GPU panel is four draws for all six fits on 16 tiles, then draws
4-31 on NGC_005_008_003 and SGC_005_014_001. A one-draw intentional stop tests
reuse without rewriting its output; a no-op resume must reuse all six fits;
an independent replay compares coarse and fine arrays exactly. One additional
real empty-retained tile is checked within the same allocation.

All checks passed:

- 720 fine draws and 240 shared coarse draws in the main panel; all saved
  checksums, dimensions, finite/positive densities, ordered eigenvalues and
  core/galaxy trace identities passed the independent output audit.
- The intentional-stop draw retained exactly the same bytes and modification
  time after resume. No-op resume reused six fits and generated zero draws.
- Independent replay reproduced all eight coarse/fine array bundles exactly,
  across both seeds (six additional fine and two coarse draws).
- The real empty tile NGC_000_009_002 passed all six fits with galaxy arrays of
  shapes (0,3) and (0,), while retaining normal field outputs (six additional
  fine and two coarse draws). There are 726 reusable fine draws in the main
  output tree, plus six replay-only draws.
- Maximum main-panel relative coarse-mass error: 7.233e-16. Maximum core trace
  absolute error: 9.936e-15. Saved float32 outputs also passed trace checks.
- Peak PyTorch allocated/reserved memory: 1.893/2.777 GiB. Slurm sampled a
  device-memory maximum of 3,382 MiB (~3.30 GiB) and GPU utilization up to 98%.
  Main-step MaxRSS was 4,023,628 KiB (~3.84 GiB). These are observed preflight
  peaks, not a guarantee for a different batching or precision configuration.

The fixed main GPU workflow took 4,394.230 s (73m14s), before its final output
audit: intentional stop 34.473 s, preview 2,287.574 s, no-op resume 12.935 s,
28-draw extension 2,008.646 s and independent replay 50.300 s. The separate
empty-tile workflow took 58.285 s. Slurm includes interpreter startup, audit,
command handoff and allocation setup/release overhead beyond these payload
timers.

The complete dense north 32-draw tile exported to a 117,422,080-byte tar archive
in 0.898 s, including verification of all 514 members. Its original bundles
remain intact. A separately verified source archive records 326 files, with
the scientific Python hashes checked against INDEX and additional launchers,
tests and configuration files recorded at closeout. This demonstrates the
archive format/export path, not sustained survey-wide I/O or archive-backed
inference resumption.

The budget calculation weights the two preflight tiles in each cap/shell by
that stratum's parent tile count. It separates coarse generation from the six
fine fits, measures write/startup overhead and reports memory. CPU preparation
is corrected for measured serialization overhead and concurrency speedup.
Final merge and Page packaging remain unbenchmarked and must retain a separate
provisional resource envelope.

Measured projections for all 2,472 tiles, including already-tested addresses:

| Work | Projection |
| --- | ---: |
| Four draws, all six fits, core computation | 97.293 GPU-h |
| Four draws, including conservative restart overhead | 99.599 GPU-h |
| Additional 28 draws, including conservative restart overhead | 697.195 GPU-h |
| Complete 32-draw campaign, core computation | 778.342 GPU-h |
| Complete 32-draw campaign, including conservative restart overhead | **796.794 GPU-h** |
| Serial survey preparation | 1.985 CPU-node-h |
| Eight-worker survey preparation, measured scaling | 0.624 CPU-node-h |
| Eight-worker preparation plus a fresh global index | 0.676 CPU-node-h |

The final conservative preview wall/compute factor is 1.02371; preparation's
wall/condition-computation factor is 1.46721. The first budget receipt used only
the resumed preview timer, producing 784.965 GPU-h. That timer omitted work
already committed by the intentional stop. `BUDGET_CLOSEOUT.json` includes the
stop's entire 34.473 s, including duplicate startup, giving the conservative
796.794 GPU-h projection above. The original receipt is preserved. The report
helper was corrected without rerunning inference; its updated source is saved
separately in `budget_closeout_source.tar`. The original source snapshot remains
the exact inference-time code record. These are engineering
extrapolations, not statistical confidence bounds. The estimates do not claim
that 16 tiles exercise all filesystem contention, extreme densities or queue
conditions. An unchanged source-bound run can reuse this index and completed
draws; no new preparation/inference for those addresses is required.

## Proposed production envelope -- not authorized or submitted

The GPU timings support retaining the earlier conservative GPU limits. CPU
preparation no longer needs the earlier two eight-hour cap jobs.

| Stage | Exact proposed allocation | Ceiling |
| --- | --- | ---: |
| Remaining preparation | One CPU node, 64 logical CPUs, `desi/regular`, 02:00:00; eight four-thread workers with disjoint explicit tile lists | 2 CPU-node-h |
| Four-draw preview | Array 0-102%8, 24 tiles/task, one A100 80 GB and 32 logical CPUs/task, `desi_g/shared`, 02:00:00 | 206 GPU-h |
| Draws 4-31 | Same 103 shards and concurrency, one A100 80 GB and 32 logical CPUs/task, `desi_g/shared`, 10:00:00 | 1,030 GPU-h |
| Final merge/QA/Page packaging | One CPU node, 64 logical CPUs, `desi/regular`, 08:00:00; provisional, not benchmarked end to end | 8 CPU-node-h |

Proposed remaining ceiling: **1,236 GPU-hours and 10 CPU-node-hours**, with
1 TB scratch headroom. A 24-tile shard projects to 0.967 h for the preview and
6.769 h for the extension; the requested 2 h/10 h limits leave headroom.
At eight concurrent GPUs, expect about 4-5 days of active inference, plus queue
time and review/packaging. Account balance and storage quota are not asserted.
All jobs would request Scratch access explicitly. GPU jobs retain
`--constraint='gpu&hbm80g' --gpus=1 --cpus-per-task=32`; no 40 GB speed or
concurrency change has been benchmarked. No automatic retries/requeues.

Before launch: obtain explicit production approval, freeze the source and
archive/staging/recovery policy, validate the production launcher/manifest
partition, and release the 28-draw stage only after full-preview coverage and
numerical checks pass. The measured budget is not permission to submit.

The un-packed full 32-draw output would contain 632,832 fine/coarse draw
bundles (about 1.27 million payload/receipt files, plus bundle directories).
Before production, freeze the archive/staging/recovery policy rather than
assuming small-file performance scales indefinitely. The archive exporter does
not delete inputs or make the inference runner read directly from archives;
those operations need an explicit production handoff and cleanup policy.

## Evidence locations

Current source-bound run:
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_preflight_20261016_v2/`.

Key receipts: `INDEX.json`, `PREPARE_*.json`, `CPU_SCALING.json`,
`STAGE_*.json`, per-draw `COMPLETE.json`, `GPU_PREFLIGHT.json`,
`MEASURED_BUDGET.json`, `BUDGET_CLOSEOUT.json`, `EMPTY_TILE_VALIDATION.json`,
`packed/NGC_005_008_003_32draws/COMPLETE.json` and
`SOURCE_SNAPSHOT/COMPLETE.json`. GPU stdout/stderr: `gpu_59834307.log` and
`empty_pack_59834307.log`.

Retained failed first attempt:
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_preflight_20261016_v1/`,
including `cpu_scaling_8/FAILED.json` and per-worker logs.

The current desktop thread's default sandbox failed before starting commands
with a PAM temporary-directory error. Reviewed native-host execution worked;
`/etc/shifter/udiRoot.conf` ownership was 0:0. No permissions/configuration,
Shifter, PAM or scheduler-review settings were changed. `graphify` is absent;
targeted source inspection was used transparently instead.
