# Loa mass–redshift completeness and historical selection

**2026-10-16 clarification:** these figures compare saved inferred-environment
catalogues, not the unrestricted Loa parent. Their low-mass deficit and
0.06-0.15 gap must not be interpreted as missing observed galaxies. See the
[parent-sample audit and new mass-selected sub-VACs](loa_mass_complete_20261016.md)
for the corrected parent distribution and controlled environment analysis.

Fresh full-footprint census, 2026-10-15. Outputs:
`/global/u2/d/dkololgi/fmsu_tokyo_20261010/loa_mass_redshift_completeness_20261015/`.
Four PNG/PDF figure pairs: common-CG15 mass–z density for full footprint and
property wedge, and two-panel historical/current selection-convention figures
for each scope. Source:
`workflows/presentation/loa_mass_redshift_completeness_20261015.py`.

## Selection and provenance

Start with all5,436,413 rows of the frozen P12-A CG15 VAC. Require SUPPORTED,
unique unambiguous CIGALE GALAXY match, finite positive mass, .15<=z<.55.
No stellar-mass, SFR, AGE, AV or FAGN cut precedes completeness estimation.
4,989,187 rows pass; every ID, RA, DEC and redshift exactly matches one row
in the registered BGS_BRIGHT targeting catalogue. Extinction-corrected
r<19.5 South /19.54 North retains4,989,181. This is the available supported,
CIGALE-matched population, not a claim that every underlying photometric
object has an environment or reliable mass.

The property wedge RA[120,160),DEC[14.5,30.6) contains221,659 rows over the
full .15−.55 range. It is broader than the .20−.30 seven-property sample.
Histograms use all selected masses; plotted range7<=logM<=13 is a display
range only. Original old-GAT CG15 matches contribute193,706 rows full sky,
7,378 in the wedge. These already passed the old FSF selection, so they
are NOT used to estimate the old mass-completeness curve.

Current curves: each dz=.01 bin, select faintest20% in r−r_limit, then
90th (and95th sensitivity) percentile of logM+.4(r−r_limit). Require100
rows/bin. Full footprint supports all40 bins. Sparse wedge endpoints may
be omitted; no extrapolation. These are empirical flux-limit diagnostics,
not guarantees against surface-brightness, spectroscopic, mass-fitting or
selection biases. The faint-end mass-to-light distribution can itself be
incompletely sampled, particularly at high z.

| Full-footprint bin | Estimated90% log mass limit | Approximate solar mass |
|---|---:|---:|
|.15−.16|10.175|1.50e10|
|.20−.21|10.508|3.22e10|
|.29−.30|10.888|7.73e10|
|.40−.41|11.231|1.70e11|
|.50−.51|11.511|3.25e11|
|.54−.55|11.648|4.45e11|

The displayed fraction above a mass limit is not the completeness percentage.
The90% refers to the limiting-mass construction, not keeping90% of observed
objects. Histograms are galaxy counts with logarithmic colour scaling.

## Historical comparison

All masses on the common main plot use CG15. Shade the old redshift interval
.01−.06; do not draw an exact historical mass boundary in CG15 units. The
actual selection used raw FSF LOGMSTAR>=9, corresponding to physical
logM>=9.305 using the historical h=.704 conversion. The companion left panel
shows that boundary and a reconstructed pre-cut main-BRIGHT completeness
curve in FSF units; its last .055−.06 bin gives9.043 full footprint and
9.000 in the wedge. The old target selection was broader than main BRIGHT;
this curve diagnoses its BRIGHT completeness rationale, not every target
programme. The old flux convention is retained and labelled in RESULTS.

The .06−.15 gap means no overlap between these saved inferred-environment
catalogues. It does NOT mean DESI observed no galaxies there. P12-A starts
at .15 because of the known low-z mock-to-matter truth-mapping failure;
this is distinct from mass completeness and not a fundamental DESI limit.
The example's BGS Faint/secondary curves are not reproduced from unsupported
samples. Our primary curve is explicitly BGS Bright.

## Interpretation

At z=.20−.30, a1e9 mass threshold cannot create a sample complete above1e9.
Many faint galaxies are absent before any analysis cut. The verified TNG
experiments show that population and smoothing changes can weaken the true
property relation, but this completeness plot alone does not establish its
amplitude or prove the missing objects are satellites. BGS Bright restricts
the observable population; it does not prohibit a strong relation in every
property, subsample or environment definition. Independent-density agreement
and same-Loa model comparisons remain the evidence against a simple model-
failure interpretation. No new inference or training; one-thread bounded
column reads and exact catalogue joins under the authorized login fallback.
