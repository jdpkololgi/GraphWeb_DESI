# Full-footprint P12-A Loa production, 2026-09-25

Status: submitted, provisional. The unresolved input-selection discrepancy is
explained in `p12a_selection_distribution_audit_20260925.md`; it is not corrected
by this run and prevents a validated science-release claim.

Root: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1`.
The128-part plan covers15,632 disjoint output cores and5,436,413 unique TARGETIDs
at0.15<=z<0.55. All6,680,148 context galaxies are available to each needed patch.
The frozen model uses20-voxel cores,48-voxel halos,alignment8 and512 untempered
joint draws per supported galaxy. Draws are per-galaxy joint eigenvalue draws,
not a joint posterior over the entire survey field.

## Scheduler and completion

- Interactive representative part0: allocation58864498.
- Array58865082: parts1–127, maximum8 concurrent single-GPU `shared` jobs,
  45minutes maximum per task,32 logical CPUs and48GB host RAM.
- Merge58865126: depends on successful array completion,1hour maximum.
- Logs: `logs/interactive_part_000.log`, `logs/part_58865082_*.out/.err`,
  `logs/merge_58865126.out/.err`.

`FULL_VAC_COMPLETE.json` is the final technical completion marker. Absence means
no completed full-survey catalogue is claimed, regardless of scheduler states.
Science release remains false even on technical success. No automatic retries.
The final FITS name is `DESI_LOA_P12A_HALO48_FULL_SURVEY_VAC.fits`; `DRAW_INDEX.json`
binds every saved per-core draw shard by hash and row count.

## Restart and identity contract

Frozen source is under `source/` with `SOURCE_MANIFEST.json`; PLAN binds every
part's input marker and global core ownership. Every worker verifies source,
model, fitted selection, schema and training-context hashes. Each core gets
seed20260925+global_core_id, independent of job order or retry.

A core writes draws and FITS through temporary files, then atomically writes its
completion marker last. Completed cores are reused only after checking input,
ownership and both output hashes. A truncated completed file is a hard error;
no silent regeneration. Incomplete temporary work may be recomputed with the
same seed. A filesystem lock prevents concurrent writers for the same part.
An exact recovery test from complete core files but missing aggregate marker
reproduced all draw and aggregate hashes. Three contract tests cover corruption,
changed input and changed ownership. The earlier canary's5615 base/response rows
match the new runner exactly;18 benchmark cores validate16542 rows.

The merge verifies every core and draw, recomputes saved means/quantiles/class
probabilities, requires each target exactly once, and compares the merged sorted
TARGETIDs to the entire input science catalogue. Unsupported rows remain present
with null posterior quantities. Completion is written only after FITS readback.

## Product flags

Bit1: unsupported/null;2: sparse z>=.45;4: boundary distance<R/h;8: distance<2R/h;
16: posterior context outside training min/max;32: provisional on every row.
The full merge additionally sets64 at z>=.35 for the unresolved selection shift.
This is an applicability warning, not a calibrated systematic-error estimate.
Intermediate per-core files have the original bits; bit64 is applied in the
full-survey merge from each row's redshift. No galaxies are dropped for this flag.

No posterior tempering, selection refit, retraining, confirmation opening or
science-release promotion is part of this execution. Upstream population/selection
closure and independent replication remain scientific follow-up work.


## First production-shard result

Part0 completed in432s:42,670 rows,42,472 supported,198 unsupported/null,
2 outside the training feature envelope (flag16). Independent validation of
every saved draw and summary passed in10s. Allocation58864498 was released.
Array58865082 is running; merge58865126 remains dependent. Full output is not
yet complete. Measured shard throughput suggests roughly2hours at eight-way
concurrency plus queue/merge time; this is an estimate, not a scheduler promise.


## Completed full-survey result

Merge58865126 completed0:0 in14m32s. FULL_VAC_COMPLETE.json now exists and
records full identity and saved-draw validation. The catalogue contains5,436,413
rows:5,404,568 supported and31,845 unsupported/null. There are125 input-envelope
flags and847,421 selection-warning flags. All rows remain provisional.
This supersedes the historical running/submitted status above.
