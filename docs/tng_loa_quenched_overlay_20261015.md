# Direct TNG–Loa quenched-fraction overlay

Plots and results: `/global/u2/d/dkololgi/fmsu_tokyo_20261010/tng_loa_quenched_overlay_20261015/`.
Primary figure `R7_tng_loa_overlay.pdf/png` compares matching R7 targets.
Supplement `R2_tng_loa_overlay.pdf/png` deliberately contrasts TNG R2 with
P12-A R7 and cannot isolate prediction error at a common scale.

All curves now use fractional equal-number quarters of tidal trace and
sSFR<1e-11/yr, not a mixture of TNG quarters and Loa cluster/void contrasts.
TNG ranks are within .25dex mass bins; Loa ranks additionally within dz=.01.
Cells need80 rows. Loa uses the established .20−.30 wedge, seven-property
eligibility, geometry mask, and explicit BRIGHT magnitude cut. Mass ceiling
is11.5 for both samples. P12-A trace is the sum of posterior-mean eigenvalues;
this is a rank statistic of a point summary, not full posterior propagation.
No new inference or retraining. TNG curves use the three-seed geometry-only
MLPs, not GAT, with frozen original predictions.

Left panels retain each population's own mass-bin distribution. Right panels
use the exact same Loa counts as reference mass-bin weights in both samples.
Common full-range support is9.5<=logM<11.5 (41,721 TNG /89,593 Loa rows);
high-mass support10.5<=logM<11.5 (9,338 TNG /68,918 Loa). The native full TNG
panel retains66,691 rows down to9.0. Support change and mass reweighting are
both included in the upper-row left/right change; it is not an additive
causal estimate for mass alone. Loa distributions are unchanged left/right
because all its eligible bins are shared with TNG.

| R7 comparison | TNG truth | Dense TNG prediction | Sparse retrained TNG | P12-A Loa |
|---|---:|---:|---:|---:|
| Full, own mass mixtures |32.73|33.75|32.64|6.71|
| Shared support, Loa mass mixture |13.28|12.93|12.76|6.71|
| High mass, own mixtures |10.10|9.70|9.68|8.30|
| High mass, same Loa mixture |7.88|7.28|7.29|8.30|

Entries: highest-minus-lowest quarter, percentage points.95% spatial intervals
for high-mass matched TNG truth6.22−9.88; Loa7.51−8.96. Full matched TNG truth
11.89−15.07; Loa6.07−7.29. Intervals exclude galaxy-formation, SED and inference
transfer uncertainty and are not independent-volume tests.

Absolute fractions remain very different: high-mass matched TNG truth ranges
85.0−92.9%, versus55.3−63.6% in Loa. Similar slopes therefore do not establish
identical populations or validate P12-A. TNG is z=0 with instantaneous SFR
and simulation stellar masses; Loa is .2−.3 with SED-derived SFR/masses.
Brightness selection, redshift evolution, satellite mix and galaxy-formation
modelling remain unmatched. Equal environment ranks do not imply equal
physical environments in the two catalogues. The point-predicted environment
itself can attenuate or reshape the relation. No model was tuned to improve
agreement with Loa properties.

Conclusion: the observed modest contrast is plausible and comparable to the
massive-TNG contrast; full-range amplitude and absolute quenched fractions
are not reproduced. Do not say P12-A is proven to recover the expected truth.

Validation:209 direct-versus-block cell-mean checks;12 independently saved
TNG curve comparisons (machine-precision agreement); exact Loa photometry ID
parity; source hashes. Curves reconstructed, not manually shifted or rescaled.
200 spatial bootstrap draws per population; TNG three-fit averaging uses
paired block draws, not averaged confidence bounds. Source script is
`workflows/presentation/tng_loa_quenched_overlay_20261015.py`.
