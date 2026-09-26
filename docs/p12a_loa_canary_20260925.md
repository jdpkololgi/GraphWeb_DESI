# Authorized Loa posterior canary and provisional VAC

User approved real DESI inference and VAC creation after conditional review and
golden replay passed. Existing Loa full data and 18 full random catalogues retain
verified source hashes. Existing graph/wedge caches have no P12 halo48 binding
and are not reused as encoder inputs. Build new count/angular/response products
from the existing catalogues; no catalogue redownload or model retraining.

Run root: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_canary_20260925_v1`.
Code: `workflows/p12a_vac/loa_trial.py`; allocation58862859, one GPU, one-hour
bound. Other E2E compute remains separate. P12 source/selection/weights frozen.

Selection: finite positive Z_not4clus, ZWARN0, DELTACHI2>=25, GALAXY from BGS
BRIGHT full HPmapcut data. Context0.1<=z<0.6 excluding0.585<=z<0.595; output
0.15<=z<0.55. Unit count weights; no automatic SYS/PIP/ZFAIL weighting. Galaxy
angular support reproduces R0; full-random GOODHARDLOC support supplies the P3b
boundary feature. Preserve cap/PHOTSYS random-domain ownership and all18 sources.
The two roles must not be conflated into a new model response channel.

Use the golden-tested17,001-point Planck18 distance lookup,5Mpc cells,40Mpc
padding,20-voxel output cores,48-voxel context and alignment8. Rebuild bounded
patches with complete neighboring context from the full selected catalogue.
No truth, predictions or DESI posterior outcomes select the canary: choose a
median-occupancy and smallest-centre-boundary-distance core with >=16 galaxies
per cap/shell,16 distinct diagnostic cores total. This is not a random sample.

Before sampling require source/checkpoint/input hashes, valid geometry, finite
fields and mock equivalence of the sparse response evaluation. Persist region
selection before model execution. Report feature min/max envelope violations
as diagnostic flags, not invented calibration guarantees or a post-hoc selection.
Unsupported rows retain TARGETID with null posterior. Save512 untempered joint
draws, quantiles and four class probabilities at eigenvalue threshold0.2. Flags
are versioned canary flags, not a silent reuse of historical production bits or
old-model prior-width thresholds. Every row is flagged provisional.

Validate identity uniqueness, finite/ordered supported draws, normalized class
probabilities and FITS readback. Assess scope/input shifts before scale-out.
The canary is a provisional research artifact, not a full-footprint VAC or
independently confirmed scientific release. Results and completion receipts follow below.


## Completed result: first real-DESI posterior catalogue

Allocation 58862859 on nid001013 completed input preparation and inference.
New context catalogue: 6,680,148 unique TARGETIDs; 5,436,413 lie in the science
redshift range. All 18 full-random maps completed. Grid padding is the frozen
40 Mpc. The geometry-selected 16 cores contain 5,615 galaxies: 5,602 supported
posteriors and 13 unsupported rows with null draws/summaries and explicit flags.
No selected posterior context lies outside the training min/max envelope; this
is a coarse support diagnostic, not evidence of matched multivariate distributions.

Product: `DESI_LOA_P12A_HALO48_CANARY_VAC.fits` under the run root above.
Each supported row has 512 ordered joint eigenvalue draws in `shards/`, means,
5/16/50/84/95 percentiles and void/sheet/filament/knot probabilities. TARGETID,
RA/DEC/Z, cap, core ownership, base eigenvalues, boundary distance and frozen
ntilde are included. All rows carry provisional bit32. Bit1: unsupported (13);
bit2: z>=0.45 (105); bit4: distance<R/h (472); bit8: distance<2R/h (1,110);
bit16: outside feature envelope (0). Bits overlap; these are versioned canary bits.

Independent readback passed: identities unique/sorted, every target appears in
exactly one draw shard, finite ordered supported draws, all unsupported posterior
entries null, and summaries/class probabilities exactly reconstructed from the
saved draws. All721 frozen source-file hashes verify. Source/checkpoint/input
hashes and small reports are mirrored in `docs/evidence/p12a_loa_canary_20260925/`.
The512-row independent mock response comparison is exact for distance/support.

The first random scan was stopped for poor strided FITS I/O; completed catalogue
files were retained. Reading contiguous whole-row blocks gives exactly the same
counts/selection in the focused test. Optimized random preparation finished in
394s; frozen inference in70s; independent product validation in2s, all exit0:0.
Non-MPI srun steps use `--mpi=none --resv-ports=0` to avoid MPI port reservations
blocking overlapping steps. These are execution fixes, not model changes.

## Next production work

V3 technical inference is complete. Before V5 launch, extend truth-free input
comparison beyond min/max (shell/cap count-to-selection ratios and conditional
feature distributions), record the intended footprint/support mask, and benchmark
representative dense and edge shards. Build restartable uniquely owned batch
shards with atomic completion markers, complete neighboring context and a merged
TARGETID census against all5,436,413 eligible catalogue rows. Do not extrapolate
full-survey runtime from these deliberately selected16 cores alone.

The current file is a bounded provisional VAC, not the full-footprint catalogue.
Independent-phase replication, observation-model robustness, class reliability
and the release/model card remain V4/V6 gates. No new confirmation data were
opened and no encoder or posterior was fitted in this run. The Loa source
catalogues are reusable; earlier graph caches were not P12-compatible inputs.

Graphify ran from cosmic_env. Compute-node Home locking failed with errno524;
a bounded login-node refresh succeeded. This was not an unavailable executable.
