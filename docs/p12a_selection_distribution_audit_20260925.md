# Loa versus P12 mock selection audit, 2026-09-25

Audited every row of Loa full HPmapcut, upstream ph006 full HPmapcut and the
exact observed ph006 catalogue used by P12. Used only licensed ph006/training
metadata; no confirmation access, changed cuts, selection refit or retraining.

## What the cuts explain

| Redshift | SGC common-footprint Loa/mock | SGC with r<19.5 and DELTACHI2>40 | NGC common-footprint Loa/mock | NGC with r<19.5 and DELTACHI2>40 |
|---|---:|---:|---:|---:|
|0.15–0.25|1.0218|1.0199|1.0222|1.0047|
|0.25–0.35|1.0758|1.0732|1.1067|1.0759|
|0.35–0.45|1.2647|1.2601|1.2672|1.2142|
|0.45–0.55|1.9423|1.9166|1.9068|1.7810|

Common footprint means the intersection of occupied nside256 pixels; this controls
the gross footprint, not all within-pixel completeness or angular selection.
The stricter cuts are diagnostic ablations, not a revised deployment definition.

- Every baseline-selected row in both catalogues passes BGS_TARGET bit1 (value2), GOODHARDLOC, GOODPRI, LOCATION_ASSIGNED and positive NOBS_G/R/Z.
- Upstream and processed mock shell counts agree exactly. All science-shell mock rows have BOX_INDEX>=0; no extra label-validity or join loss accounts for the deficit.
- DELTACHI2>40 removes 1,853/137,328 (1.35%) of the highest-shell Loa baseline. Mock columns have no DELTACHI2 or SPECTYPE; simulation galaxies do not have a spectral-classification model.
- All selected mock R_MAG_APP values are <19.5. Loa has northern-photometry targets extending to r=19.54. Tightening Loa to19.5 removes5,078/95,587 (5.31%) in the highest NGC shell, almost none in SGC. This is a photometric-convention difference, not demonstrated erroneous target membership.
- Even the aggressive diagnostic MASKBITS==0 changes counts by only a few percent and does not close the discrepancy. It is not proposed as a replacement for the registered HPmap cuts.
- Loa WEIGHT_ZFAIL means are approximately1.001–1.004 across the eight bins. Adding mock spectroscopic failures would reduce mock counts further, not supply the missing high-z galaxies. These weights are diagnostics, not inserted into frozen count fields.

## Fibre assignment and sentinel redshifts

Full mock tables have sentinel redshifts for unassigned targets. Counting failures
using their Z_not4clus or TRUEZ would silently exclude them. The independent
MOCK_ASSIGNMENT_AUDIT instead joins full-table TARGETIDs to parent observed/RSD
redshifts before binning. Assigned and ZWARN0 counts are identical; assignment
retention is approximately77.7–84.1%. In the highest shell the full target tables
contain26,219 SGC and59,277 NGC galaxies, compared with
41,741 and95,587 successfully selected Loa galaxies. Even assigning every mock
target cannot bridge the excess. This bound concerns these supplied HPmapcut
target catalogues; it does not establish equivalence of upstream targeting.

The earlier SELECTION_AUDIT_V2 TRUEZ counters describe only available valid
full-table redshifts and must not be interpreted as total assignment efficiency.
Use the parent-joined assignment report for that question.

## Interpretation and production disposition

No accidental extra cut was found in our mock-processing chain. The tested
quality/magnitude/footprint/assignment choices do not explain the high-z excess.
Existing ph002–006 highest-shell totals range69,293–71,331 versus Loa137,328;
ph006 is not an outlying phase. The remaining explanation is unresolved upstream
mock population/photometric evolution/targeting or observation-model mismatch.
Do not claim an HOD error, missing galaxies, or a proven causal explanation yet.

Full-footprint inference may produce a provisional model-conditional research
VAC under the user’s existing authorization. It cannot be promoted to science
release from these diagnostics. Keep the frozen selection unchanged and add
QUALITY bit64 for the unresolved z>=0.35 selection discrepancy (32 remains set
for every provisional row). This flags known applicability concerns; it does not
marginalize their uncertainty or certify the lower-redshift domain.

Follow-up before science release: inspect parent/mock photometric and n(z)
construction, compare matching targeting/photometric conventions and angular
response, then test any scientifically justified observation-model correction
on mocks. Changing the fitted selection/input definition requires matched
encoder-summary/posterior validation, not a silent DESI-only rescaling.

## Execution evidence

Input distribution and selection reports live alongside the canary artifacts.
The full production ownership plan contains15,632 cores,128 balanced shards,
and5,436,413 unique science-range galaxies with complete neighboring context.
The18-core dense/edge benchmark has16,542 rows. Its saved draws validate, all
5,615 earlier-canary base predictions and response values match exactly, and
restarting after core completion reproduces the catalogue and draw hashes.
Three restart-contract tests pass. Production submission is recorded separately.

## Submitted full-footprint provisional build

Array58865082 covers parts1–127, at most8 concurrent single-GPU shared jobs;
part0 is the representative frozen-source run in allocation58864498.
Dependent merge58865126 waits for array success, then validates all128 parts
and the exact full TARGETID census before writing FULL_VAC_COMPLETE.json.
Each worker has45minutes maximum; merge has1hour. No automatic retries.
Source and runtime contracts are frozen at
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1/`.
Submission is not completion; no full-survey product is claimed yet.
