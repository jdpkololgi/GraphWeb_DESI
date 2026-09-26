# Full-survey DESI posterior diagnostics

Production and merge completed successfully. The final catalogue has5,436,413
rows, with5,404,568 supported posteriors and31,845 null/unsupported rows.
The plotting run verified the full FITS SHA against FULL_VAC_COMPLETE.json,
checked finite normalized probabilities and nonnegative interval widths, and
verified hashes for every individual-example draw shard. No fitting or new
posterior sampling was performed. CPU allocation58870808 was used for plotting.

[Four-page PDF](figures/p12a_loa_full_20260925/DESI_posterior_diagnostics.pdf).
Reproduction: `workflows/p12a_vac/plot_posteriors.py ROOT NEW_OUTPUT_DIRECTORY`.
Source hash, example TARGETIDs and numerical arrays: accompanying PLOT_SUMMARY.json.

## What is plotted

1. Survey overview: distributions of per-galaxy posterior medians (not pooled
   posterior densities), central68% widths against redshift, mean class probabilities
   and distributions of maximum class probability. Histogram display clips only
   the combined0.1–99.9% median range; denominators retain every galaxy in each group.
   Width bands describe16–84% population spread, not uncertainty on a mean.
2. Equal-area sky pixels at0.15<=z<0.25: mean filament probability and mean λ2
   interval width; nside64 pixels require at least10 supported galaxies. Radial
   projection is not a3D field reconstruction. Coordinate labels are ICRS.
3. Four galaxies' marginal posteriors from their512 saved joint draws. One per
   shell selected nearest the median λ2 width among interior supported objects
   without envelope flag16. These examples are deliberately selected, not random.
4. Pairwise joint draws for the lowest/highest-redshift examples, showing
   dependence between ordered eigenvalues. No independent-marginal multiplication.

## Numerical reading

Median maximum class probability is0.74414;40.01% of supported galaxies have
maximum probability>=0.8. These probabilities have not been calibrated against
DESI environmental ground truth. Mean probabilities over the observed supported
galaxies are20.45% void,41.43% sheet,30.63% filament,7.49% knot. They are not
volume fractions, nor a selection-corrected population estimate.

Median68% interval widths are0.1587,0.1747,0.2119 for λ1/λ2/λ3. In the first
z bin[.15,.175), they are.1399/.1510/.1846; in the last[.525,.55),
.4801/.5643/.6561. The increased width is observed, not proof that all systematics
are covered. Mean filament/knot probabilities rise at high z while void/sheet
probabilities fall. The differing tracer sample, weaker constraints and known
selection mismatch prevent interpreting that trend as cosmological evolution.

All figures remain provisional. No DESI coverage/SBC/TARP can be obtained from
these plots without environmental truth; mock coverage is separate evidence.
The z>=.35 applicability warning is visible, and release qualification remains
open. All four exported PNGs were visually inspected for readable labels and
clipping; the PDF contains the same four rendered figures.

![Overview](figures/p12a_loa_full_20260925/01_posterior_overview.png)
![Sky projection](figures/p12a_loa_full_20260925/02_sky_projection.png)
![Individual posteriors](figures/p12a_loa_full_20260925/03_individual_posteriors.png)
![Joint posteriors](figures/p12a_loa_full_20260925/04_joint_posteriors.png)
