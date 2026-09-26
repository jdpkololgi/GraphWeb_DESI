# Selection closure: historical targeting and Kibo/Loa — 2026-09-26

Diagnostic-only work authorized by the user. No catalogue selection changes,
repairs, retraining, inference or VAC replacement. Allocation58885256,
cosmic_env. The work narrows the cause but does not establish the exact raw
mock generator recipe or full real/mock completeness equivalence.

## Kibo and Loa are both DR2-era reductions

[DESI DR2 Results II, sectionII](https://arxiv.org/html/2503.14738v2#S2)
describes Kibo as the homogeneous reduction of observations through9April2024,
and Loa as its rerun after fixing a co-addition error. Roughly0.1% of measured
redshifts changed significantly across the overall dataset. DR1 used Iron,
not Kibo. These public aggregate numbers are not substituted for our BGS audit.

Our model currently uses the kibo-v1 ph006 mock. A Loa mock also exists:
`/global/cfs/cdirs/desi/survey/catalogs/DA2/mocks/SecondGenMocks/AbacusSummitBGS_v2/altmtl6/loa-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.h5`.
Its HDF5 compression requires hdf5plugin registration; the initial read failed,
then succeeded after importing the already-installed package. No input modified.

Full real Kibo-v1/v1 and Loa-v1/v2.1 catalogues, plus the Kibo/Loa ph006 mocks,
were counted on the same fixed common sky used in the previous count plot.
Comparisons mix LSS catalogue versions as well as spectroscopic productions;
they are not a clean measurement of the co-addition change alone.

| z | Real Loa/Kibo SGC/NGC | Mock Loa/Kibo SGC/NGC | Real Loa / Loa mock SGC/NGC |
|---|---|---|---|
| .15–.25 | 1.003 / 1.005 | 1.045 / 1.037 | .977 / .985 |
| .25–.35 | 1.002 / 1.005 | 1.045 / 1.038 | 1.029 / 1.066 |
| .35–.45 | 1.003 / 1.004 | 1.045 / 1.037 | 1.210 / 1.221 |
| .45–.55 | 1.002 / 1.004 | 1.046 / 1.036 | **1.856 / 1.839** |

The historical Kibo/mock mismatch is therefore not DR1 versus DR2. Matching
Loa products reduces the count discrepancy modestly but preserves its steep
redshift dependence. No decision to substitute Loa mocks has been made.

The tile unions encoded in TILES contain5164 for real Kibo and mock Kibo,
5165 for real Loa,5171 for mock Loa. Real Loa adds25529; mock Loa adds22245,
22855,22856,25257,23823,24563,25529. These are tile unions represented in these
catalogues, not universal lists of every usable exposure. Historical mask and
bad-exposure equivalence is not implied by equal tile counts.

Real successful-ID intersection:5,391,797; Kibo-only23,271, Loa-only44,616;
62 shared successful IDs change redshift by>.001. Mock successful intersection:
4,967,614; Kibo-only594, Loa-only198,232; zero such redshift changes among shared
IDs. Larger mock count changes reflect catalogue/retention differences, not
changed redshifts for shared objects. Evidence: PRODUCT_CROSSWALK.json.

## Closure point1: selected-DESI targeting provenance and replay

All5,436,413 science TARGETIDs join to the original
`target/catalogs/dr9/1.1.1/targets/main/resolve/bright` catalogue. All source
file headers identify desitarget1.1.1. A historical real fibreassign header
(tile23286) independently points to that directory. Extracted actual function
bodies from the desitarget1.1.1 Git tag, including isBGS, notinBGS_mask,
isBGS_colors and isBGS_sga. Used its verified BRIGHT/CLUSTER geometry bits1/13.

**Every joined galaxy passes historical BRIGHT selection; zero target-bit
mismatches.** Gaia/REF_CAT fields are now included;5166 have SGA-like REF_CAT.
None fail the Gaia test. Stored original g/r/z fluxes equal Loa fluxes exactly
for every science galaxy, so updated fluxes do not explain this discrepancy.
This closes the missing-fields/revision check for the selected DESI sample.
It does not establish completeness of the parent galaxy population, which
includes objects never targeted or without successful spectra.

Evidence: TARGET_REPLAY.json, historical cuts and geometry source snapshots.
The replay is read-only; its purpose is verification, not reselecting the VAC.

## Closure point2: mock production lineage — partial

Located an actual archived ph006 per-tile fa shell script and output headers.
The script for tile20090 uses the historical real sky/secondary/footprint files,
observation rundate2021-05-14 and specified assignment margins. However it names
fiberassign4.0.0 whereas the saved fba file records5.7.2.dev3588 and desitarget
2.8.0.dev5597. This may reflect a rerun, wrapper or changed environment: the
script alone cannot certify the execution. No causal effect is inferred.
Headers archived in MOCK_HEADERS.json; source script remains at
`altmtl6/Univ000/fa/MAIN/20210514/fa-020090.sh` under the mock root.

forFA6, pota-BRIGHT and combined-assignment headers have no generator command,
Git revision or LF/K-correction parameters. Raw CutSky header contains area
metadata but no recipe. The nearby docs inventory describes archival storage,
not generation. The current local preparation code and published BGS method
remain supporting context, not verified production receipts for these files.
Exact HPmap input hashes, bad-fibre/exposure lists and full catalogue command
versions remain unpinned. We have not silently treated matching filenames or
current LSS defaults as historical proof.

## Closure point3: retention and conditional population comparisons

Rejoined every Kibo/Loa full mock TARGETID to the official ph006 forFA BRIGHT
parent and recovered RSDZ even for unassigned rows. Every BRIGHT parent has
R_MAG_APP<19.5. This avoids treating sentinel unassigned redshifts as real z.

High-shell(.45–.55), fixed common sky:

| Stage | SGC | NGC |
|---|---:|---:|
| Raw CutSky r<19.5, previous census | 30,124 | 65,661 |
| forFA BRIGHT parent | 29,801 | 65,218 |
| Kibo full HPmap targets | 26,197 | 59,251 |
| Kibo assigned / ZWARN0 | 21,246 | 49,826 |
| Loa full HPmap targets | 26,636 | 59,887 |
| Loa assigned / ZWARN0 | 22,216 | 51,637 |
| Real Loa successful galaxies | **41,243** | **94,985** |

All mock full targets pass GOODHARDLOC/GOODPRI; assignment and ZWARN0 counts
coincide. Kibo assignment retention in the high shell is81.1%/84.1%, Loa83.4%/
86.2%. Neither later angular cuts nor assignment can explain why even the
upstream raw bright population is already smaller than observed DESI successes.
This argument is conditional on the still-unproven equivalence of apparent
magnitude definitions and parent photometry.

Existing conditional diagnostics remain relevant: uniform numericr<19.5,
cap separation, and matching redshift composition inΔz=.01 retain a high-z
colour difference≈.17/.16mag between DESI and mock. Galactic cap is not PHOTSYS:
SGC is southern photometry here; NGC mixes both systems. Real fibre-magnitude
and Gaia selection were now replayed. Mock forFA has no equivalent noisy
fibre/Gaia/SGA observables, so a direct identical-observable completeness replay
cannot be manufactured. Real unobserved targets lack measured redshifts;
parent completeness as a function of true z is not identifiable by simply
binning the failed/unobserved real catalogue's reported redshift.

## How photometry and population modelling can make a redshift trend

Write apparent magnitude schematically as
`m_r = M_r(z) + DM(z) + K_r(z, SED)` with consistent units/conventions.
Then the flux limit implies `M_limit = m_limit − DM − K`. A too-positive mock
K-correction or too-faint mock luminosity evolution excludes galaxies that
would otherwise pass. At high z the cut probes the steep bright tail of the
luminosity distribution, where a small magnitude error removes a much larger
fraction. A passband mismatch can itself depend on z and galaxy SED, rather
than being a constant zeropoint offset. Colour modelling affects both that
mapping and which galaxy types cross the flux limit. This is a selection effect;
colour is not directly fed into our current position/count encoder.

Quantitative sensitivity from the existing raw-parent census: moving only the
diagnostic magnitude threshold from19.5 to19.6 raises counts by9.6–9.7% in the
first shell, but47.0–47.8% in the last. A.04mag change raises them by3.8% versus
17.3–17.4%. These are finite-difference illustrations, not proposed repairs or
estimates of an actual offset. A deficient evolving luminosity function can
also reduce high-z counts without any implementation error. The observed mock
colour offset alone does not determine the sign or amplitude of the r-band
K-correction error; the actual recipe must be recovered.

## Current disposition

Point1 selected-sample replay is closed. Points2/3 progressed substantially but
remain incomplete for the actual generator recipe and full population-dependent
selection model. Next evidence required is the official raw-CutSky generation
receipt/code/config and executed full-catalogue mask/observation recipe. This
is a provenance requirement, not permission to apply a correction. No repairs
were made. Keep the existing VAC provisional.

Artifacts: `docs/evidence/p12a_closure_20260926/`; scripts
closure_catalogue_crosswalk.py, closure_targeting_replay.py,
closure_mock_retention.py, closure_summary.py. Independent stage counts match
crosswalk counts exactly, all TARGETID joins/replay checks pass. All three
successful analyses completed on allocation58885256; initial HDF5 read failure
is retained in its log and was resolved by loading the required installed filter.
