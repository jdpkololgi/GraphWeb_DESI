# ph006 LSS selection test — 2026-09-26

User now authorizes creating a test phase and iterating generation/selection
parameters toward Loa agreement. This supersedes the previous diagnostic-only
restriction for this isolated experiment; no production VAC replacement or
training is included. ph006 is already exposed. No reserved phases are read.

1. Run the pinned official LSS Y3 BGS preparation on existing ph006 CutSky,
   with seed62026, rbandcut19.5, downsampling n, apply_mask n, isProduction n.
   Compare generated BRIGHT rows/IDs/photometry to official forFA6_nomask.fits.
   This tests preparation reproducibility; current code is not automatically
   the historical code. FAINT RNG and TARGETID ordering need explicit checks.
2. Build an isolated common-sky development table from the existing raw phase,
   retaining raw-row and halo linkage. Compare Loa quality25 and official-style
   quality40, including numerical bright limits and PHOTSYS. Never give noiseless
   mocks invented DELTACHI2, SPECTYPE or fibre photometry.
3. Sweep a bounded two-parameter apparent-magnitude mapping, separately from
   the unresolved physical meaning of stored M and LF evolution. Fit only SGC
   count shape, check NGC and colour/magnitude distributions. Coarse grid then
   one bounded refinement, all trials retained. No arbitrary n(z) resampling.
4. Existing Loa mock/raw shell retention may be used only as a clearly labelled
   fixed-retention screening approximation. It is not a rerun of fibre assignment
   or a validated quality model for newly promoted galaxies. Count-matching
   candidates must subsequently regenerate targeting/assignment/observation and
   pass joint population/clustering checks before being accepted.
5. If no candidate passes, report the mismatch and missing modelling inputs;
   do not tune endlessly or silently flatten real evolution. Matching Loa used
   in tuning is calibration, not independent validation.

Working root: /pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_mock_test_ph006_20260926_v1.
Allocation58895657, CPU/cosmic_env, one hour. Config is
configs/p12a_mock_test/ph006_20260926.json. Large trial catalogues stay on Scratch;
code, compact evidence, logs and decisions are committed here and in SCIENCE_LOG.

## Completed first iteration

The official upstream preparation completed in 4m29s after setting DESIMODEL
explicitly to the installed main data tree. The initial attempt and missing-data
error are retained; no shared environment or CFS catalogue was changed.

### Baseline reproduction and footprint discrepancy

- Official archived forFA: 18,062,024 total, 10,540,172 BRIGHT.
- Current pinned Y3 preparation: 18,419,903 total, 10,747,480 BRIGHT.
- Exact RA/DEC/RSDZ join recovers all 10,540,172 archived BRIGHT galaxies; no
  archived galaxy is missing. Shared TRUEZ, R_MAG_APP, R_MAG_ABS, G_R_OBS and
  BGS_TARGET are identical. TARGETID shifts follow extra rows before ID assignment.
- All 207,308 extras are outside IN_Y5; archived BRIGHT has zero outside IN_Y5.
  Thus the current source's bypassed Y5 cut is a concrete reproduction difference.
- On our fixed common sky the extras are only 77 SGC / 141 NGC galaxies in
  0.45 <= z < 0.55: 0.258% / 0.216% of archived forFA BRIGHT counts. This
  difference cannot explain the roughly 84–86% Loa/mock high-shell discrepancy.
- This content comparison establishes preparation behaviour, not the historical
  Git commit, FAINT random state, executed fibre-assignment version or all masks.
  LSS preparation preserves the raw magnitude columns for shared galaxies; it
  does not introduce the unresolved -.8(z-.1) relation.

A sky-only join was deliberately stopped at repeated directions. The successful
join uses unique RA/DEC/RSDZ triples; no many-to-many or first-match approximation.

### Parameter screening

Prepared 16,637,115 raw common-sky science-range galaxies and 5,309,765 successful
Loa galaxies with numerical 12 <= r < 19.5. Baseline raw and full-Loa shell counts
reproduce the previous independent census exactly. Compared 30 coarse and 19
additional refined trials of m_trial = m_raw + delta_m - extra_q*(z-.1).
This is an observer-magnitude sensitivity model; extra_q is not an identified
physical LF evolution parameter and does not redefine stored absolute magnitude.

Best broad-shell SGC count screen: delta_m=+0.075 mag, extra_q=0.6. It uses fixed
baseline Loa mock/raw retention in dz=.01 bins. Ratios of predicted counts to
Loa successful quality25 counts are:

| z | SGC | NGC |
|---|---:|---:|
| .15–.25 | 1.011 | 1.018 |
| .25–.35 | 1.041 | 1.032 |
| .35–.45 | 1.052 | 1.085 |
| .45–.55 | .976 | 1.041 |

These broad shells conceal residual fine-bin shape. No tested candidate is
within 10% in every dz=.01 bin of both caps. Even selecting the existing grid
post hoc by fine-bin SGC error leaves extrema of 26.6% SGC / 16.7% NGC.
The broad-shell candidate's high-z numerical g-r discrepancy worsens from
mock-minus-Loa -0.163 to -0.191 mag in SGC, and -0.152 to -0.180 in NGC.
Passbands/flux estimators remain unaligned; these are not estimates of an
intrinsic physical colour error. It would be improper to tune intrinsic colour
until the observer-band mapping is established.

The Loa quality40/quality25 count ratio is only 1.005–1.013 across broad shells.
These are different combined rules: quality25 includes SPECTYPE=GALAXY;
quality40 uses DELTACHI2>40 without that condition, so it can retain more rows.
Neither change explains the large redshift trend. No mock DELTACHI2, SPECTYPE,
fibre flux or noisy spectral-success field was invented.

A 7,352,166-row trial parent FITS is saved on Scratch with raw-row/halo linkage,
original magnitudes and separate R_MAG_TRIAL/BGS_TARGET_TRIAL columns. It is
explicitly TESTONLY, has no fresh fibre assignment and is not a fitted-quality
mock. Existing production inputs and the VAC remain unchanged.

### Disposition and next iteration

The first magnitude-only parameter family is not accepted. Do not propagate
its count fit into training or the VAC. Continue from the reproducible baseline:

1. Pin the Y5 toggle in the test recipe and quantify remaining targeting/HPmap
   masks and assignment retention separately from photometric generation.
2. Resolve SDSS-Petrosian versus Legacy observer photometry and the raw stored-M
   evolution definition, then vary physically coherent LF/colour/SED parameters
   through that observer model. The failed numerical colour comparison is a
   diagnostic, not a licence to fit an arbitrary colour offset.
3. Establish the spectral-success observation model for relevant Loa cuts;
   rerun assignment after any target-population change. Fixed retention is not
   adequate to certify a changed targeting catalogue.
4. Require fine-z joint magnitude/colour/size/sky agreement and clustering, then
   truth-known conditional coverage on separately allocated confirmation data.
   SGC tuning and the already-exposed NGC check do not provide that confirmation.

Evidence is in docs/evidence/p12a_mock_test_20260926; the standalone figure is
in docs/figures/p12a_mock_test_20260926. The 16 lightweight GraphWeb P12 tests pass.

Full Y5 replay closes the BRIGHT preparation population check: applying IN_Y5
and the documented sequential postcut TARGETID assignment reproduces every
checked column in the same order for all10,540,172 BRIGHT objects. This is
stronger than equal counts or approximate positional matching. See Y5_REPLAY.json.
The seeded test baseline is materialized separately under y5_baseline; its FAINT
random realization and later observation pipeline are not historically reproduced.

Materialized Y5 test baseline:18,062,788 total targets,10,540,172 BRIGHT. The
764-row difference from archived total is in the new seeded FAINT realization.
Independent read-back checks10,000 spaced BRIGHT rows against the archived file
exactly across IDs, positions, redshifts and magnitudes, plus row count/last ID.
An initial last-row read used unsupported negative fitsio indexing; rerun with
explicit positive row index passed. No catalogue change was needed.
