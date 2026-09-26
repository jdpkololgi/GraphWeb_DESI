# Selection diagnostics protocol — 2026-09-25

Authorized by the user's count-plot and controlled-test request. Allocation
58871404. Frozen halo48 encoder/posterior; no fitting or VAC replacement.

1. DESI and ph002–006 counts in 40 bins over 0.15–0.55 on their common
   NSIDE256 occupied-sky intersection, split by cap. Compare frozen expected
   geometric counts; five-phase range is not a systematic uncertainty interval.
2. Reuse the existing 50k ph006 evaluation rows and 512 draws with natural
   weights. Compare identical-row truth classes, encoder classes, posterior
   probabilities, and fine-redshift coverage. ph006 is exposed development data.
3. Reuse the eight previously chosen golden cores. Global Bernoulli thinning
   probability 0.5, RNG seed20260925; choose up to64 retained, supported galaxies
   per core with the same RNG. Evaluate those same galaxies in every arm:
   baseline; unchanged counts with half expected density; half counts with
   original expected density; half counts with half expected density.
   Expected-density change applies both to the field contrast and posterior
   log-ntilde covariate. Keep the angular response/window fixed to isolate
   density effects. Query galaxies are not artificially protected from thinning.
   Require rebuilt baseline fields and predictions to match the stored export
   at atol=rtol=1e-5 before interpreting perturbations. Use512 posterior draws
   and paired RNG seeds per core/arm. Report paired probabilities, eigenvalue
   shifts, interval coverage and Brier score against unchanged known truth.
   These dense, geometry-selected cores are not representative coverage tests.
4. Trace upstream catalogue/photometry selection without changing cuts after
   seeing results. Preserve reserved phases and existing production outputs.

No flat-redshift environment-fraction requirement; flux selection can alter
population mix. Matching n(z) alone does not qualify simulation-to-DESI transfer.
