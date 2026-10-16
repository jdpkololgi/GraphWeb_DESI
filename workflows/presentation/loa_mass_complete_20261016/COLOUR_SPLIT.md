# Observed-colour split

`colour_split.py` extends the completed parent audit without changing its outputs.
Red means strictly `0.75 < GR_OBS < 1.5`; blue means strictly `GR_OBS < 0.75`.
Exact boundaries and colours >=1.5 are excluded from both plotted groups;
nonfinite colours are excluded. GR_OBS uses Milky-Way-corrected observed fluxes,
with no rest-frame or K correction. Colour is not a quenching label.

Reconstruct colour for the registered Bright parent with exact TARGETID joins,
then apply the same faintest-20% / 95th-percentile limiting-mass estimator
separately to red and blue populations. Each redshift bin retains the maximum
of the previous all/quiescent/starforming limit and both new colour limits.
A fixed threshold is the maximum across the selected redshift interval,
rounded upward to 0.01 dex. This remains empirical flux completeness.

Compare zmax=.30,.25,.23,.22 with zmin=.20. Before inspecting property outcomes,
select the smallest tested zmax below .30 having at least 20,000 galaxies in
each colour and at least 100 in each hard environment class within each colour.
Retain .20-.30 as a baseline. This count rule is a pragmatic support threshold,
not an uncertainty or power guarantee. Actual matched support is reported.

For each property, require common .1 dex mass / .01 redshift strata across all
red/blue and soft/hard environment curves, with Kish count >=20 per curve.
Repeat at .05 dex. The stellar-mass panel controls redshift only. Each property
has its own valid rows. Shared reference weights permit comparisons among all
eight curves on their retained support; unsupported strata are not imputed.
No star-forming or atomic-gas sequence is assumed or fitted.

Use 256 paired HEALPix nside=8 bootstrap draws and show conditional 95% intervals.
These exclude SED errors, completeness-boundary uncertainty, survey-selection
systematics and simulation-to-data transfer uncertainty. Red/blue membership
may reflect redshift, dust and stellar populations, not only star formation.

Output FITS preserve all cached input columns and add COLOUR_CLASS (-1 excluded,
0 blue, 1 red). Existing SELECT_FIXED020030 and LOGM_LIMIT95 retain their original
meaning; use the new file's MCFIX/MCZLO/MCZHI header for its stricter colour-audited
selection. Large products are written to a new Scratch directory.

Validation independently reconstructs percentile order statistics and every
reported standardized mean, and roundtrips every source column at up to 257
rows per FITS. A manifest records code hashes, parent manifest hash and outputs.

## Observed support and follow-up

No tested window passed the predeclared 20,000-blue / 100-per-hard-class gate.
The original .20-.30 sample has only 375 blue galaxies; .20-.22 has 1,922.
The latter is retained as an explicitly exploratory lower-mass demonstration,
not promoted to a gate-passing primary result. `colour_followup.py` reports raw
curves and controls within each colour separately, preserving the failed joint
analysis in v1. Soft-only controls and joint soft/hard controls are distinguished;
both .1 and .05 dex strata are retained. Within-colour controls do not enforce
the same mass/redshift reference between red and blue.

The .20-.30 blue sample has no controlled support. At .20-.22, soft-only controls
retain 1,462 blue galaxies at .1 dex, 996 at .05 dex; joint soft/hard controls
retain only 447 in one .1 dex stratum and none at .05 dex. Red controls retain
almost the whole population. Sparse binary outcomes can have degenerate
nonparametric bootstrap intervals when a class has no observed events; a zero
interval in that situation is not evidence of a population probability of zero.

`colour_export.py` exports final presentation figures and all raw/controlled
means and cluster-minus-void contrasts to CSV. Original v1/v2 manifests remain
unchanged; a separate presentation manifest binds the final rendered figures.

## Requested 0.15–0.20 extension

`colour_lowz.py` and `colour_lowz_export.py` prepare the same colour audit and
within-colour analysis for 0.15<=z<0.20. The input is the existing variable-cut
sub-VAC: the new fixed threshold must exceed every old bin limit in this range,
so it cannot omit any qualifying object through the input selection. Fresh
original-VAC parity is checked in 257 rows through SOURCE_ROW. Outputs use a
new `loa_colour_lowz_20261016_v1` directory. Execution completed after explicit CPU approval on allocation 59857020.
The fixed limit is 10.60, with 518,962 red and 14,909 blue galaxies. The
20,000-blue gate still fails. Results, sensitivity checks and validation are
recorded in `docs/loa_colour_lowz_20261016.md`; the allocation is released.
