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
