# Population-model tests: ANY selection, P/Q, versions and Uchuu

**Superseded for execution ordering by [Mock-to-LOA Alignment Investigation](plan_mock_to_loa_alignment_investigation_20260926.md).** Retained as the earlier design record.

User follow-up 2026-09-26. This protocol expands the active VAC closure work.
Header inventory completed; the experiments below are not yet executed. Production
VAC and model weights remain unchanged. Use exposed phases for development;
freeze choices before accessing any reserved confirmation phase.

## Newly verified candidates

Raw Abacus ph000 at
`/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/`:
v0 has63,696,561 rows; v0.1 has63,897,781; v1 has74,748,522.
v1/z0.200 contains only ph000 in the current inventory. Its schema replaces
HALO_INDEX with HALO_ID, includes box/replicated coordinates and velocities, and
lacks the stored IN_Y masks. More total rows does not establish more high-z
BRIGHT galaxies or a matched footprint. Inspect recipe/header provenance, rebuild
the same angular selection, and compare counts/photometry before any inference.
ph000 is training-exposed, not independent confirmation.
Header evidence: evidence/p12a_reference_mocks_20260926/ALTERNATIVE_INVENTORY.json.

Public Uchuu BGS exists and is readable at
`/global/cfs/cdirs/desi/public/edr/vac/edr/uchuu/v1.0/BGS-BRIGHT_Uchuu/`.
Inspected S64:79,952 rows, Z, apparent/absolute magnitudes, rest colour, K and E
corrections, NZ and FKP; no delivered halo IDs. This is a viable population check,
not yet a compatible environmental training set. The official documentation
identifies102 SV3 footprints and BGS0<z<0.6:
https://data.desi.lbl.gov/doc/releases/edr/vac/uchuu/ .

DR1/iron and EDR/SV3/fuji data are alternatives, not additional independent Loa
observations. DR1 is a subset of the main-survey observing interval; restrict Loa
to common TARGETIDs/footprint for a reduction check. SV3 has a different observing
strategy and small volume. Public availability is documented at
https://data.desi.lbl.gov/doc/releases/ . Local Y1/iron and SV3/fuji LSS paths
are present. DA2 Kibo/Loa and Loa LSS v1.1/v2.1 are processing alternatives,
not interchangeable new independent surveys.

## A. BGS_ANY candidate for our full science sample

1. Use the deepest available raw parent, retain original photometry and halo IDs,
   and establish sufficient parent density in every z/sky bin. A threshold cannot
   supply galaxies absent from the parent; raw r20.2 truncation can limit changes.
2. Define target expected counts for exactly our Loa footprint and quality sample.
   Fit a smooth M_lim(z), with selection M<M_lim(z), using the generator's explicit
   absolute-M and evolution convention. Include expected assignment/redshift
   success when solving for the parent threshold. Do not fit each observed n(z)
   fluctuation exactly: that can absorb actual large-scale structure.
3. Tune on SGC, evaluate NGC without retuning; both are development comparisons,
   not independent mock confirmation. Record finite-volume uncertainty and cap
   photometric differences. Keep a flexible ANY threshold as a diagnostic control
   alongside a physically constrained LF/photometry model.
4. Rerun preparation/assignment for changed targets with the intended BRIGHT
   observation priorities; do not inherit observed FAINT assignments or multiply
   by old BRIGHT retention. Recompute response/expected-density fields.
5. Count fitting is an admissible empirical calibration, not intrinsically
   improper. Qualification additionally needs conditional magnitudes/colours,
   luminosity-dependent clustering, host/satellite mix and posterior calibration.
   A threshold-only candidate can be useful even if stored mock apparent
   photometry is imperfect, but must be labelled as an empirical sample mapping.

## B. Identify the physical mapping error and test P/Q directly

First reproduce the existing catalogue magnitude mapping and cuts with its
actual recipe. Recover the pivot redshift, luminosity-function parameters,
cosmology/h units, K-correction bands and whether M already contains evolution.
Avoid applying the evolution term twice.

On identical exposed-phase halos, geometry and stochastic seeds, run:
- baseline P1.8/Q0.7;
- P-only changes at fixed Q and K/photometry;
- Q-only changes at fixed P and K/photometry;
- K/observer-photometry changes at fixed P/Q;
- then a joint P/Q model, constrained by DESI luminosity-function evidence.

Vary the actual target LF/remapping, not just a plotting E-correction, and do not
approximate P with arbitrary weights on the final observed points. Recover the
baseline exactly before interpreting response derivatives. Screen modest
parameter changes, then bound the range using LF constraints rather than
selecting unrestricted parameters only from n(z). If a deeper halo/galaxy parent
is required, distinguish parent regeneration from remapping existing objects.

Record n(z), luminosity functions, r and g-r quantiles at fixed z/PHOTSYS, fibre
magnitude where physically defined, and selected host/satellite distributions.
Compare same-galaxy SDSS/Legacy photometry against colour and z; use compatible
FastSpecFit K-corrections and completeness-aware LF estimation. If P/Q improves
counts but worsens colours or luminosity-dependent clustering, reject the claim
that it alone solves the mismatch. A response to P/Q establishes sensitivity;
causal attribution needs baseline provenance and successful joint predictions.

## C. Compare genuine alternatives before dismissing them

First screen raw Abacus v1 versus v0.1 on exposed ph000 with matching geometry
and selection; version labels alone do not prove a physical improvement.

Run the public Uchuu ensemble against its intended SV3 BGS sample with identical
magnitude/redshift cuts, comparing completeness-corrected data to intrinsic mocks
or explicitly simulating observation losses. Do not compare intrinsic Uchuu raw
counts directly with incomplete Loa successes. Use the ensemble to quantify
small-footprint scatter. Subsequently test a verified full-sky/DR2 parent under
Loa observation selection; existing accessible SV3 products cannot cover Loa by
simple rotation. Continue locating the actual DR2 BGS parent and provenance.

Use DR1 versus Loa common-area/common-selection comparisons to separate release
processing from changing population coverage. Earlier successful DR1 BAO tests
used a different absolute-magnitude threshold; do not infer that full-BRIGHT
DR1 automatically solves our high-z population mismatch.

## Gates and priority

Priority: raw-v1 screen and matched Uchuu/SV3 population test; exact baseline
recipe replay; factorial P/Q/K tests and smooth ANY candidate; then fresh
assignment, joint clustering/population checks and posterior validation.
No demand for exact n(z) equality in independent realizations. Predefine
uncertainty-based tolerances before viewing confirmation results. Existing
confirmation reservations remain in force. Do not migrate training, shrink the
science redshift range, or replace the VAC on count agreement alone.
