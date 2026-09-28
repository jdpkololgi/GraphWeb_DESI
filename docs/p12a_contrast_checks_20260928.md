# Frozen-model environment–property contrast checks

2026-09-28. The existing P12-A VAC, selection, mocks and posterior values are
unchanged. These checks find informative mock predictions, appreciable contrast
dilution, and a measurable CIGALE association. They do **not** establish a large
improvement in galaxy-property separation over the saved June model or validate
real-DESI coverage.

## Evidence

- [Three-page figure PDF](figures/p12a_contrast_checks_20260928/contrast_checks.pdf)
- [Mock discrimination](figures/p12a_contrast_checks_20260928/01_mock_discrimination.png)
- [Same-galaxy model comparison](figures/p12a_contrast_checks_20260928/02_same_sample_models.png)
- [Neighbour-density reference](figures/p12a_contrast_checks_20260928/03_density_reference.png)
- [Mock numbers](figures/p12a_contrast_checks_20260928/MOCK.json),
  [observational numbers](figures/p12a_contrast_checks_20260928/OBSERVATIONS.json),
  [run](figures/p12a_contrast_checks_20260928/RUN.json),
  [artifact hashes](figures/p12a_contrast_checks_20260928/ARTIFACT_SHA256.json).
- Reproducer: `workflows/p12a_vac/contrast_checks.py`; CPU allocation 59022935.

## 1. Existing-mock contrast recovery

Only the already-exposed ph006 50,000-row evaluation and its saved 512 posterior
draws per galaxy were read. No new inference, fitting, or confirmation-phase
access occurred. Dataset/checkpoint provenance and the softplus ordered-eigenvalue
inverse transform were checked; source hashes are saved. Classes count physical
eigenvalues above 0.2; natural weights restore the sampled mock population.

For the 21,052 evaluation galaxies at 0.2 <= z < 0.3, argmax accuracy is 76.45%.
The multiclass Brier score is 0.327 versus 0.681 for a spatially cross-fitted
constant-prevalence baseline (lower is better). All four class-specific Brier
skills are positive. However, true-knot recall is only about 46%: about 48% of
true knots are assigned to filaments. Global accuracy hides this weakness.

An explicitly synthetic property `y = true_class / 3` has a unit true
knot-minus-void contrast. Probability-weighted classifications recover 0.680
(spatial-bootstrap 16–84% interval 0.674–0.688); argmax recovers 0.839
(0.833–0.848). Thus class mixing can suppress this imposed contrast by about
32% or 16%, respectively. This is a diagnostic of attenuation, **not simulated
SFR**, and is not a correction factor for the DESI measurements. One independent
Gaussian-property realization provides a null diagnostic, not a null ensemble.
The ph006 bootstrap does not measure variation between independent mock phases.

## 2. Identical CIGALE galaxies: current versus June predictions

The saved CIGALE wedge is joined by unique TARGETID to the current VAC with
coordinate checks. The eligible wedge contains 94,415 supported galaxies with
positive finite mass and SFR, 120 <= RA <= 160, 14.5 <= DEC <= 30.6 and
0.2 <= z < 0.3. No FastSpecFit properties are used. The historical comparator is
specifically the cached June 2026 `desi_wedge_flowjax_linear_si` posterior at
threshold 0.2, whose historical join averaged duplicate graph-node predictions.
It does not represent every earlier transductive model.

All four estimators below use the **same 61,739 galaxies in 68 common cells**,
with widths 0.1 dex in log stellar mass and 0.01 in redshift. Each cell must have
effective membership >=20 in every class for every estimator. Cell means use
the same pooled galaxy-count reference distribution; soft weights are posterior
membership probabilities. Low sSFR means log10(sSFR/yr^-1) < -11.

| Estimator | Knot minus void, mean log sSFR (dex) | Low-sSFR fraction difference |
|---|---:|---:|
| Current probabilities | -0.238 | +9.07 percentage points |
| Current argmax | -0.302 | +11.63 percentage points |
| June probabilities | -0.247 | +9.36 percentage points |
| June argmax | -0.277 | +10.60 percentage points |

The current probability-weighted contrast has interval -0.249 to -0.221 dex;
the June counterpart -0.258 to -0.229 dex. Paired differences and all intervals
are in the JSON. There is no dramatic gain over June in this property diagnostic.
Similar bimodal distributions can coexist with a roughly nine-percentage-point
difference in low-sSFR fraction; visual overlap alone misses that association.

Intervals use paired resampling of 22 HEALPix nside=8 sky blocks, 128 attempted
draws (123 valid for the fine comparison). Draws with empty retained class/cells
are excluded. These are 16–84% spatial-bootstrap intervals, not full uncertainty
including CIGALE systematics, mass/SFR errors, selection, or model transfer.
Control bins and the pooled reference distribution remain fixed in resampling.

Requiring current max probability >=0.8 leaves 44,777 initial galaxies but only
2,692 in five cells shared by all four estimators. Larger contrasts there are
not representative evidence of improvement: the usable population is different.
Standard coarser controls are also recorded. A single within-fine-cell
permutation changes the soft low-sSFR contrast from +7.77 to -0.12 percentage
points on its own 89,264-row common support. This supports an association beyond
those finite control bins, but is not a permutation-derived significance test.

## 3. Same-galaxy reference independent of neural weights

Count neighbours excluding self in a redshift-space sphere of radius
7/0.6766 = 10.34585 Mpc, using Planck18 comoving coordinates and all relevant
full-VAC tracers, including those outside the analysis wedge. The padded region
contains 316,480 tracers. Normalize counts by sphere volume times the frozen
local `NTILDE_MPC3`; require boundary distance greater than the sphere radius.
This leaves 91,839 interior queries. Local expected density approximates its
integral over the sphere. The estimator shares the survey and selection
calibration with the model, but uses no neural predictions.

Density quartiles are constructed within dz=0.01, without artificially splitting
ties. They are **not** void/sheet/filament/knot truth, and top-hat tracer density
is not the Gaussian-smoothed matter tidal field. No independent geometric
web-finder such as DisPerSE was run.

On 65,566 identical galaxies in 77 fine control cells:

| Estimator | Endpoint mean log sSFR contrast | Low-sSFR fraction contrast |
|---|---:|---:|
| Current probabilities, knot minus void | -0.237 dex | +9.18 pp |
| Current argmax, knot minus void | -0.298 dex | +11.68 pp |
| Density Q4 minus Q1 | -0.247 dex | +9.29 pp |

The density contrast interval is -0.257 to -0.236 dex and +8.75 to +9.76 pp.
Expected model class and normalized neighbour density have Spearman rho=0.787
on the interior sample. The reference does not reveal an enormous property
contrast hidden solely by the model. Different class definitions and windows
prevent interpreting the endpoint comparison as an accuracy ranking.

## Interpretation and integrity

The measured trend is not absent. Mock class mixing attenuates contrasts, and
knots deserve particular caution. The density reference and historical comparison
do not support claiming a dramatic new improvement in property separation.
None of these tests identifies physical truth for individual DESI galaxies.
All panels use internally matched samples; samples differ between panels.

The VAC SHA256 before and after is
`cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc`.
The joined diagnostic cache and its hash are recorded in OBSERVATIONS.json.
Unit checks cover common-reference composition removal, paired identical
estimators, injected within-cell contrast recovery, and the eigenvalue transform.
The provisional-release limitations remain: real-DESI coverage is unverified,
training-population mismatch is known, especially above z~0.35 and severely at
0.45<z<0.55, and redshift trends cannot be interpreted as physical evolution.
The lower-redshift wedge is not automatically certified by these checks.
