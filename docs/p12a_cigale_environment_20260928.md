# Current Loa environments and saved CIGALE wedge properties

Follow-up: [cosmic-web literature assessment](p12a_property_literature_assessment_20260928.md)
compares these effects with other finders and distinguishes physical weakness
from possible inference attenuation. Visual similarity alone is not a model-quality verdict.

The current P12-A VAC shows a moderate, monotonic association with CIGALE galaxy
properties in the old sky wedge. The association is weaker after controlling
redshift and stellar mass, but does not disappear. This is an observational
consistency check, not a measurement of true-class population relations or
proof of calibrated DESI environment probabilities.

Per the user's clarification, this analysis uses **CIGALE masses and SFRs only**;
no FastSpecFit property results are presented. The attempted full-survey join
was abandoned. Neither the original selection, saved posteriors, trained model
nor mocks were changed. No mock-repair work was restarted.

## Figures

[Four-page PDF](figures/p12a_cigale_properties_20260928/LOA_CIGALE_environment.pdf)

1. [Property means and low-sSFR fractions](figures/p12a_cigale_properties_20260928/01_cigale_environment_trends.png): raw versus standardized relations.
2. [Relations within stellar-mass bins](figures/p12a_cigale_properties_20260928/02_fixed_mass_relations.png): redshift-balanced means in 0.5-dex mass bins. Some residual mass variation within a bin remains possible.
3. [Mass–sSFR distributions](figures/p12a_cigale_properties_20260928/03_mass_ssfr_distributions.png): probability-weighted selected-sample distributions with common axes/colour scale. These panels are not mass/redshift standardized.
4. [Membership sensitivity](figures/p12a_cigale_properties_20260928/04_membership_sensitivity.png): probability weights, argmax labels, argmax with max P>=0.8, and CIGALE CG5.

## Primary results

All primary results weight each galaxy by its **current VAC** environment
probability. Low-sSFR means log10(sSFR/yr^-1)<-11; this is a descriptive
quiescence proxy, not a directly verified quenched label.

| Property | Void | Sheet | Filament | Knot | Knot minus void |
|---|---:|---:|---:|---:|---:|
| Mean log10 stellar mass, common redshift mix | 10.698 | 10.743 | 10.784 | 10.812 | +0.113 dex |
| Mean log10 sSFR, common redshift and mass mix | -10.645 | -10.711 | -10.783 | -10.850 | -0.205 dex |
| Low-sSFR fraction, common redshift and mass mix | 43.48% | 45.88% | 48.63% | 51.25% | +7.77 percentage points |

The spatial-bootstrap 16–84% intervals on the knot-minus-void contrasts are
[0.108, 0.117] dex for mass, [-0.217, -0.194] dex for log sSFR, and
[7.33, 8.35] percentage points for the low-sSFR fraction. These intervals do
not include property-model errors, transfer systematics or calibrated joint
uncertainty in the environment field. They are not full error budgets.

Before adjustment, the knot-minus-void log-sSFR contrast is -0.469 dex and
low-sSFR fractions rise from 39.55% to 55.96% (+16.40 points). Thus differing
mass/redshift populations explain a substantial part of the raw association.
The controlled -0.205 dex is a factor 0.624 in the **geometric mean** sSFR,
not a ratio of arithmetic mean SFRs. It is a meaningful but not enormous effect.
Mass is the stronger driver apparent in the fixed-mass panels. No matched
comparison with the months-old model was made, so these figures do not establish
that P12-A improves on that model.

## Sample, join and provenance

The cached property source is
`/pscratch/sd/d/dkololgi/graphweb_desi/flowjax_inference_outputs/desi_wedge_cigale_hz/desi_wedge_env_props.parquet`.
It contains the historical CIGALE CG15/CG5 masses and SFRs. Its original upstream
`/global/cfs/cdirs/desi/users/manasvee/1_prepped_data/DESI/loa/` location is now
permission-denied; this run pins the surviving cache, not a newly verified
upstream catalogue. The cache was described in SCIENCE_LOG on 2026-07-09.

Only TARGETID, positions, redshift, MASS_CG/SFR_CG and MASS_CG5/SFR_CG5 were read.
**No historical environment probabilities, eigenvalues or hard labels were
read or reused.** Current probabilities came from the existing full-survey
P12-A halo48 FITS by unique exact TARGETID joins.

Of 111,503 cache rows, 103,814 match the current VAC. Three disagree in redshift
by |delta z|/(1+z)>=0.001 and are excluded; all position checks pass within
1 arcsec. Six further matched rows fall outside the current-coordinate wedge.
These analysis exclusions, including exact IDs and discrepant redshifts, are
listed in [RESULTS.json](figures/p12a_cigale_properties_20260928/RESULTS.json).
The intended display sample is RA 120–160 degrees, Dec 14.5–30.6 degrees,
0.20<=z<0.30, with supported current posteriors and finite positive CIGALE
mass and SFR. This gives **94,415** galaxies. Nonpositive SFR is excluded,
not automatically interpreted as quenched. No missing-value imputation was used.

Valid-property availability among supported matched wedge galaxies is
90.83%, 90.99%, 91.10%, 90.99% by probability-weighted void/sheet/filament/knot.
This near-equality does not establish that the goodPhoto/property-selected
sample is unbiased or complete. The unmatched original cache rows and the rest
of the Loa footprint are outside the scope of this comparison.

The original VAC SHA256 was checked before and after analysis and remains
`cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc`.

## Estimators and limits

For raw class k means use sum_i(P_ik y_i)/sum_i(P_ik). These are descriptive
soft-membership summaries. They are not automatically unbiased estimates of
E[y|true class k], even if individual class probabilities were calibrated.

Standardization uses 0.025 redshift bins and, for sSFR/fraction, 0.25-dex CIGALE
mass bins over 8<=log10 M*<13. Retain only cells where every environment has
at least 20 effective objects, N_eff=(sum w)^2/sum(w^2). The reference cell
mixture is the pooled retained sample and is identical for all four classes.
The adjusted sSFR/fraction sample retains **94,028** galaxies in 32 cells;
mass needs only the four redshift cells and retains all 94,415. Raw and
standardized curves consequently differ slightly in support as well as weighting.
Finite bin widths leave some residual confounding possible.

Uncertainty uses 96 bootstrap resamples of the 22 occupied HEALPix nside=8
sky regions, holding the common-cell reference mixture fixed. These are
approximately 7-degree angular blocks, not independent mock realizations.
Finite sky area, unequal regions and correlations across boundaries limit the
interpretation. Galaxy-property measurement uncertainties and transfer errors
are not propagated. No p-value or coverage claim is made.

Sensitivity panels require common support separately for each estimator; their
retained samples differ. Stronger contrasts among high-confidence labels can
reflect sample selection as well as reduced membership dilution, and do not
make that subset the preferred population. CG5 uses CG5 mass in its controls.

All inherited VAC limitations remain: real-DESI coverage unverified; known
training-mock mismatch especially above z~0.35 and severe at 0.45<z<0.55;
redshift trends must not be interpreted as physical evolution without independent
validation; this lower-z wedge is not automatically certified. This analysis
supplies a property cross-check only, not science-release qualification.

## Reproduction

Source: `workflows/p12a_vac/property_environment.py`. Use cosmic_env:

```bash
python workflows/p12a_vac/property_environment.py --test
python workflows/p12a_vac/property_environment.py \
  --out docs/figures/p12a_cigale_properties_20260928
```

Use an absent output directory. Run the catalogue operation on a CPU allocation.
Validated with allocation59021332/nid004160. Synthetic checks verify constant
properties, removal of a pure population-composition effect, and preservation
of a known within-stratum difference. Four PNG figures and the PDF were visually
reviewed. Source, cache, VAC and output hashes accompany the numerical results.
