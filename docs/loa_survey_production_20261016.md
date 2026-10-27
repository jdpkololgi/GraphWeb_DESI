# Six-fit Loa survey campaign: staging, recovery and launch contract

Status: **four-draw GPU preview array59850809 submitted**,103 tasks with throttle8;
currently queued. Dependent one-hour CPU audit59850820 is pending
`afterok:59850809`. **The32-draw extension is not submitted.** CPU preparation
allocation59849839 on `nid200002` completed0:0 in19m23s and is released.
All2,472 tiles are prepared;4,628,192 rows are retained and4,238,281 meet the
defined property eligibility mask. All28 focused tests pass on compute, the
real32-draw archive round-trip passes514 member hashes/no-op modification times,
and the real24-fine/eight-coarse-draw audit passes on a recovered tile.

The initial CPU allocation request was rejected by auto-review; the user then
explicitly approved the exact request. No scheduler workaround or security-setting
change was used. Preview and dependent audit submissions were accepted normally.

## Authorized scope and resource accounting

The user requested the full six-fit run within1,236 GPU-hours and10 CPU-node-hours,
at most eight concurrent GPUs, with a four-draw preview reviewed before extending
to32. This is inference from existing checkpoints, not retraining, mock
calibration, a new physical operator, or Page publication.

| Stage | Planned maximum resource reservation | Status |
| --- | --- | --- |
| Lifecycle validation and preparation | One CPU node,64 logical CPUs,desi/interactive,2h |59849839 completed0:0,19m23s,released |
| Four-draw preview |103 single-GPU jobs x2h,desi_g/shared,array0-102%8 |59850809 submitted/queued |
| Preview audit/review |One CPU node,64 logical CPUs,desi/regular,1h |59850820 submitted; afterok59850809 |
| Additional28 draws |103 single-GPU jobs x10h,desi_g/shared,array0-102%8 | Not submitted; preview-review gate |
| Final merge/packaging |One CPU node,desi,7h | Provisional reserve; not submitted |

GPU jobs request one80GB A100 and32 logical CPUs, Scratch license, no automatic
requeue. Total maximum GPU reservation is206+1,030=1,236 GPU-hours. CPU reservations
sum2+1+7=10 node-hours. Walltime limits are ceilings, not forecasts or allocation
balances. Actual usage must be recorded from Slurm accounting. Failed attempts
count against the same budget; no automatic resubmission is licensed.

Preparation used0.32306 CPU-node-hours. Slurm accounts256 allocated logical CPUs
for this exclusive CPU node; the payload step requested64. GPU submission
readback verifies `ArrayTaskThrottle=8`, `gres/gpu=1`,32CPUs,2h, `Requeue=0`,
account `desi_g`, feature `gpu&a100&hbm80g`; requested `shared` routes to
`gpu_shared`. The CPU allocation is absent from the active queue. Submission
timestamps in Slurm's display are local Pacific time, not UTC. Queue start time
is unknown; the earlier submission dry-run was only an estimate, not a reservation.

The preflight projection was99.6 GPU-hours for four draws and697.2 more for32,
including conservative restart overhead. It did not benchmark whole-survey merge.
At quota inspection on16October, Scratch used10.86/20.00TiB and1.03/10.00M inodes;
this was a point-in-time check, not reserved capacity. Shared single-GPU jobs and
resource charges follow the [NERSC QOS policy](https://docs.nersc.gov/jobs/policy/).

## Fixed scientific inputs

Use the2,472-tile preflight INDEX and its checksum-bound catalogue without
repeating the expensive cross-match. Retained interpolation-compatible rows:
4,628,192 of5,436,413;808,221(14.87%) are excluded at tile interpolation edges,
before support/galaxy-property cuts. All tiles, including49 with zero retained
galaxies, remain in the inference plan. These cuts must be visible in downstream
coverage reports, not presented as full parent-catalogue completeness.

Six fits are legacy CFM, expanded U-Net endpoint and transformer velocity,
each seeds17/29; within seed, coarse draws are shared. NFE128 and addressed
tile/seed/draw noise remain unchanged. Unit-density exterior beyond the finite
local-plus-wide domain and independently generated tiles remain explicit
limitations. More tiles/draws do not establish calibrated or coherent survey
posteriors. No mock truth, oracle coarse field or sealed phase is used.

New campaign root only:

    /pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1

Original preflight, earlier pilot, catalogue inputs and checkpoints are read-only
sources. Nothing in those directories is cleaned up. No files are transferred
to third-party services; archives contain private galaxy data and are not Pages.

## Archive and staging contract

1. **Prepare once.** Eight workers own disjoint tile lists and write atomic
   prepared bundles. Keep the prepared arrays and receipts for all tiles.
2. **Fixed partition.** INDEX order is partitioned into103 disjoint shards of24.
   A shared/exclusive index-lock contract and per-tile/per-shard locks prevent
   concurrent writes to the same output. At most eight shards execute at once.
3. **Atomic draw checkpoint.** Existing unchanged inference writes each coarse
   or fine bundle via temporary directory then rename. Receipt binds source,
   inputs, checkpoint, tile/noise address and payload checksum. Partial hidden
   directories are not successful draws and are not automatically deleted.
4. **Commit archive.** After all six fits finish a tile, pack its prepared data,
   coarse draws and fine draws in one uncompressed tar (NPZ members are already
   compressed). Four draws gives66 members;32 gives514. Verify exact membership,
   every member hash, receipt/payload linkage, tile/seed/draw identities,
   checkpoint/coarse linkage and INDEX identity. Flush payload/receipt/directory
   before atomic publication; preserve container checksum and manifest.
5. **Remove redundant staging only.** Restrict deletion to the literal new
   campaign root with its matching CAMPAIGN marker. First verify the complete
   archive and every existing file in the proposed removal list. Unlink only
   expected checksum-matched coarse/fine arrays and receipts, then remove empty
   draw/model directories. No recursive directory deletion; no prepared files,
   archives, unknown files or hidden incomplete directories are deleted.
   Record each cleanup and its recovery archive. Original preflight is untouched.
6. **Preview remains recoverable.** Retain the four-draw archives even after a
   separate32-draw archive is made. There is no automatic archive garbage
   collection, quota cleanup, purge workaround or HPSS/CFS write.

Only active or failed shards retain raw draws; successful shards promptly
reduce to prepared inputs plus archives. Eight simultaneous24-tile shards
bound normal active raw staging to192 tiles, instead of all2,472. Failed partial
shards remain for diagnosis and need explicit retry/budget decisions.

## Recovery and failures

- Verify container, exact member allowlist, hashes and provenance **before**
  restore. Never call unrestricted tar extraction; reject links/path escapes.
- Restore missing files through flushed temporary files and no-clobber links.
  Verify existing files and preserve their modification times on no-op recovery.
  A mismatch fails closed; no original or existing conflicting file is overwritten.
- The next inference checks the restored bundles with its original bindings.
  Draw count, worker and shard do not change the noise address.
- Interrupted cleanup is recoverable: the archive is already complete and
  missing redundant files can be restored. Interrupted restore fills only
  missing files; existing mismatches block. Partial hidden files remain diagnostic.
- A complete archive takes precedence for a completed tile; a later approved
  retry can finish redundant cleanup without regenerating that tile.
- Source/checkpoint drift, bad hashes, incomplete preparation, numerical
  failures and expired walltime stop the affected job. Do not silently skip
  problematic tiles, rewrite receipts, resample favorable outputs or retry.
- Scratch archives are durable relative to job termination, **not backups**:
  Scratch purge/loss risks still apply. Long-term archival is a separate action.

The implementation uses unchanged live scientific source with startup/repeated
hash checks and saves exact source tar snapshots, including untracked lifecycle
files. Those sources must remain unchanged while the campaign runs. No Git
commit or push is implied by a source snapshot. A changed implementation needs
an explicit new binding/version, not editing the active campaign in place.

Git revisions at preparation launch are GraphWeb_DESI
`68cd7010a3ff73cf0313b24a8ce622edc86597a8` and Illustris
`f5d16cd0617249e132d11cbb6c1d1573b29ebdb1`. Untracked implementation files are
included in the explicit source snapshots, so these revisions alone do not
identify the execution source. The326-file scientific/source snapshot tar has
SHA256 `76ecacf8c87fce6f530b46857f31adf5fd8622b53003e0001ec15b446f9d9c0e`;
the seven-file production snapshot tar has
`2e01832a135b580cd319cd4b39e44b2f247bb42a7b212950ba21242ccaee32e6`.

## Validation and preview review

Ten standard-library fixture tests pass in0.249s on login: exact member counts,
round-trip/no-op/cleanup/recovery, wrong-root refusal, validate-before-delete,
interrupted cleanup, preservation of conflicting files, corrupted archive,
path traversal/symlink rejection, wrong membership and INDEX mismatch.
Shell syntax checks pass. On CPU59849839 the same ten plus nine runner tests
pass in0.796s. Nine additional preview-auditor tests pass in0.038s. The prior
preflight's inference/replay and empty-tile numerical checks remain separately
recorded.

Real32-draw archive recovery passed on CPU59849839: all514 members matched and
a second restore preserved modification times. Receipt `RECOVERY_VALIDATED.json`
binds source archive SHA256
`7d2f398932167aa6748b836b4cc1e8561e03ed2b65246ae2968d76d1f8418e95`.
All2,472 tiles finished preparation with source-bound completion receipts before
GPU launch. `PREPARATION_COMPLETE.json` records the row counts above and INDEX
SHA256 `5b30ad1261b7da27a8ad6065aabf8576c73abc794b065b6d2c3bc617269a375a`.
The real validation output is retained; the original archive is never altered.

`PREVIEW_AUDITOR_VALIDATED.json` verifies24 real fine draws, eight coarse draws,
saved-core block mass against the decoded coarse field, array shape/finite/
positive/order/trace checks and class-count reconciliation on recovered tile
NGC_005_008_003 (10,407 retained /10,269 eligible rows). All source catalogue
columns on the fresh prepared version of that tile are exactly equal to their
catalogue slices. The NumPy warning concerns only
`CIGALE_SPECTYPE`'s HDF5 ASCII encoding annotation; it does not alter the checked
string values. This is a representative value check, not an independent reread
of every catalogue row. The real audit took4.10s including normalizer/checkpoint
loading, and the four audit-source files are saved with tar SHA256
`295e31baa67dcf9fc0e977ed5f64a54300524ef084b5c963ef015e6693332129`.

After the full four-draw preview, review all103 shard receipts and2,472 archives,
59,328 fine draws and19,776 shared coarse draws. Check missing/duplicate tile
ownership, hashes, source identity, numerical invariants, array dimensions,
empty/sparse cases, cap/redshift/boundary coverage, eligible-row counts, extreme
field summaries, throughput, storage and actual budget consumption. A clean
scheduler exit alone is insufficient. Observed-data review cannot test posterior
calibration and must not reclassify the four draws as a calibrated estimate.

Before extension, additionally require all preview tasks to be terminal in
Slurm accounting, reconcile their actual GPU-hours and any failures against the
ceiling, and confirm no preview GPUs remain allocated. The throttle of eight is
per array; never overlap the preview and extension arrays and thereby exceed
the campaign-wide limit. A payload receipt alone does not prove job release.

The extension script refuses to run without `PREVIEW_REVIEW.json` explicitly
authorizing extension and binding the current campaign manifest. This receipt
is not created automatically by success of either the preview array or its CPU
audit. Only after review
may draws4…31 be generated, reusing exact draws0…3. The extension remains
unsubmitted at this document's current status. No background agent wakeup has
been configured; the deterministic dependent audit can complete without a chat,
but interpretation, budget reconciliation and extension submission require a
subsequent review in this thread. If any preview task fails, the afterok audit
does not run; no failed-task retry or automatic extension is enabled.

## Repository entrypoints

- `workflows/presentation/loa_survey_production_20261016.py`: initialization,
  source checks, preparation, archive verification, restore, scoped cleanup, shards.
- `workflows/presentation/test_loa_survey_production_20261016.py`: lifecycle tests.
- `workflows/presentation/pack_loa_survey_tile_20261016.py`: archive publication.
- `workflows/presentation/loa_survey_preview_20261016.slurm`: reviewed preview array.
- `workflows/presentation/loa_survey_extension_20261016.slurm`: gated extension array.
- `workflows/presentation/audit_loa_survey_preview_20261016.py`: complete preview
  archive/array/ownership/eligibility audit; never creates extension authorization.
- `workflows/presentation/loa_survey_preview_audit_20261016.slurm`: optional
  afterok-dependent CPU audit, one hour from the reserved review budget.
- `docs/e2e_fresh_eyes_briefing_20261016.md`: scientific methodology/draw-count handoff.
