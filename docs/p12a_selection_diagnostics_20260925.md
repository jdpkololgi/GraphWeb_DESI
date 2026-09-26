# DESI/Abacus selection and posterior diagnostics — 2026-09-25

The count mismatch is real and increases with redshift. It is not a survey-wide
factor of two. The filament/knot redshift trend, however, also exists in mock
truth; it must not be attributed wholly to selection mismatch or called physical
evolution. All production VAC rows remain provisional.

## Matched-sky counts

[Count plot](figures/p12a_selection_20260925/01_matched_counts.png),
[PDF](figures/p12a_selection_20260925/01_matched_counts.pdf),
[data](figures/p12a_selection_20260925/COUNTS.json).

DESI and ph002–006 use exactly the same NSIDE256 occupied-sky intersection:
4,177 deg² SGC and9,544 deg² NGC (rounded). Forty bins, Δz=.01. Unit galaxy
counts, successful-redshift selections pinned to existing production. No
truth-label requirement for this count comparison. Five phases share a mock
population prescription; their min/max range is not a systematic-error interval.

| Redshift | DESI / mean Abacus, SGC | DESI / mean Abacus, NGC |
|---|---:|---:|
| .15–.25 | 1.014 | 1.032 |
| .25–.35 | 1.074 | 1.096 |
| .35–.45 | 1.256 | 1.265 |
| .45–.55 | 1.973 | 1.911 |

Pooled common-footprint ratio:1.086812 (+8.68%). Last .54–.55 bin ratios:
4.757SGC and3.992NGC; the coarse highest-shell factor≈2 conceals a steep tail.
The excess exceeds the five-phase spread in the high-z bins, but no independent-
galaxy Poisson significance is claimed. Frozen expected counts integrate
ñ(z)r²dr times common solid angle: geometric approximation, not the exact
apodized voxel normalization. Its highest-shell DESI ratios are1.952/1.871.
The expected curve follows mocks closely until the last few bins, where smoothing
softens their sharper decline. The mismatch is therefore not just a fitted-n(z)
normalization mistake.

## Mock truth versus posterior trends

[Comparison](figures/p12a_selection_20260925/02_truth_encoder_posterior.png),
[coverage](figures/p12a_selection_20260925/03_mock_coverage.png),
[data](figures/p12a_selection_20260925/PH006.json).

Reused50,000 saved ph006 evaluation rows and512 draws per row. Apply saved
natural weights to truth, encoder hard classes and mean posterior probabilities
on the same rows. ph006 is already exposed development data, not confirmation.
Class threshold is0.2, ordered eigenvalues, fixed target epochz=.2.

At z=.15–.175, ph006 truth filament/knot fractions are.253/.054; at .525–.55,
.464/.310. Corresponding mock posterior means are.260/.058 and.455/.292.
DESI posterior means are.280/.068 and.477/.292. Thus the strong redshift trend
already exists in the mock-selected galaxy population. The encoder's hard
classes alone diverge substantially at high z; posterior probabilities repair
much of this aggregate difference. A mean posterior probability is not an
unbiased population-PDF estimate, nor is matching it a DESI calibration test:
the posterior shares the mock prior and redshift/selection conditioning.

Last bin has only168 evaluation galaxies; central68% empirical coverage is
.714/.649/.637, versus.700/.717/.711 in the first bin. These fine bins are
exploratory and spatially correlated; they do not replace registered coverage
gates. No real-DESI truth coverage has been measured.

## Controlled density perturbations

[Protocol](p12a_selection_test_protocol_20260925.md),
[plot](figures/p12a_selection_20260925/04_selection_sensitivity.png),
[summary](figures/p12a_selection_20260925/SENSITIVITY.json).

All eight previous golden cores passed baseline field and stored-encoder parity
at atol=rtol=1e-5.468 identical retained, supported query galaxies,512 draws
per arm, common RNG, frozen angular response. Global Bernoulli thinning avoids
artificially retaining query galaxies. No targets, weights or posterior refit.

| Arm | Mean filament/knot P | λ2/λ3 central68% coverage | Brier |
|---|---|---|---:|
| Baseline | .479/.170 | .669/.639 | .413 |
| Same counts, half expected density | .555/.230 | .470/.365 | .567 |
| Half counts, original expected density | .411/.113 | .485/.521 | .482 |
| Half counts, half expected density | .489/.174 | .650/.654 | .446 |

Underestimating expected density moves λ2/λ3 medians by about.50/.68 of their
baseline68% interval width and worsens calibration on these paired galaxies.
Consistent density adjustment restores much of the aggregate baseline behaviour
for random thinning, with residual information loss. The altered ñ at fixed z
is an intentional stress test, potentially outside the training joint feature
support; it does not estimate DESI bias or demonstrate that empirical DESI n(z)
can simply be substituted. Eight dense cores are not a representative calibration
sample. Raw outputs and per-core metrics:
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_selection_sensitivity_20260925_v1`.

## Upstream trace

[Upstream plot](figures/p12a_selection_20260925/05_upstream_counts.png),
[census](figures/p12a_selection_20260925/UPSTREAM.json).

Streamed all64,058,709 rows of the registered ph006 raw CutSky, restricted to
the identical common sky. Raw r<19.5 counts at z=.45–.55 are30,124SGC and
65,661NGC, before any fibre assignment. DESI successful observations on that
sky number41,243/94,985: even this raw-parent upper bound is deficient by
factors1.369/1.447. This localizes a substantial discrepancy upstream of our
observed-truth join, fibre assignment and post-redshift cuts.

Using Z_COSMO instead of RSD Z changes these raw high-shell counts by−0.56%/
+0.62%. Applying IN_Y5 barely changes them. There is no hard raw redshift cutoff
at.55: the census continues through.60 and source rows extend to≈.8. The raw
magnitude limit is20.2, fainter than the BRIGHT cut. Our parent builder selects
R_MAG_APP<19.5 then matches forFA identities; it has no .55 redshift cut.
The observed-truth builder selects finite positive Z_not4clus and ZWARN0;
its selected counts already match the upstream LSS table exactly.

Raising the diagnostic raw limit from19.5 to19.54 increases high-shell counts
by≈17%, but still leaves them below DESI successes even before assignment.
Larger changes can strongly alter the tail; these are not validated photometric
corrections. The current records do not establish equivalence between mock
R_MAG_APP and Loa extinction-corrected FLUX_R photometry, nor identify the
upstream luminosity-evolution/K-correction prescription as the specific cause.
No HOD error is proven. Prior full-table cut audit showed stricter DELTACHI2,
magnitude and hardware cuts do not explain the discrepancy.

## Disposition and next work

1. Preserve the completed VAC and its high-z/provisional flags. Do not flatten
   environment fractions or silently substitute a DESI n(z) fit.
2. Obtain/pin the actual generator recipe behind these raw apparent magnitudes:
   luminosity function, evolution, K-correction and photometric-system mapping.
   Compare DESI/mock colour and magnitude distributions conditional on z and
   cap; distinguish tracer-population changes from random thinning. The source
   census narrows the cause but does not close this physical selection question.
3. Register a physically justified matched-population/observation candidate.
   Rebuild fields and selection consistently, then representative mock tests of
   coverage, conditional bias and posterior predictive observables. Use exposed
   development phases for choice; preserve independent phase reservations.
4. If the changed fields move encoder summaries, regenerate OOF summaries and
   posterior fits, with encoder retraining contingent on measured performance.
   Validate on an eligible independent mock before producing a versioned VAC
   replacement. Matching n(z) alone is insufficient.
5. Science release also needs truth-free DESI closure against density/clustering
   and external observables with selection controls. Mock coverage remains
   conditional on the mock model; exact DESI coverage is not directly observable.

Execution: allocation58871404, cosmic_env. Counts, saved-draw analysis,32
core/arm runs, raw-parent census and plots completed successfully. No reserved
phase opened, no fitting, no production artifact altered. Source metadata and
frozen manifest hashes are recorded; source FITS hashes were not recomputed in
this diagnostic run. Figure inspection and count/identity checks completed.
