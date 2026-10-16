# Observed colour distributions after joint mass/redshift matching

All original CIGALE mass selections and observed Milky-Way-corrected Legacy
Survey g−r colours are retained. Environment predictions are fixed P12-A
P_WEB. No FastSpecFit masses or rest colours enter this analysis.

[Primary highest-probability-class figure](figures/loa_matched_colour_20261016/matched_colour_primary_hard.png)
shows unweighted distributions on common support (top), matched distributions
(middle), and differences from the equal-environment mean matched distribution
(bottom). Each has a PDF counterpart. Alternatives:
[primary probability weights](figures/loa_matched_colour_20261016/matched_colour_primary_soft.png),
[finer hard matching](figures/loa_matched_colour_20261016/matched_colour_fine_hard.png),
[finer probability weights](figures/loa_matched_colour_20261016/matched_colour_fine_soft.png).

## Estimand and matching

Primary cells are 0.10 dex in log stellar mass × 0.01 in redshift. Retain only
cells with effective sample size ≥20 in every environment under BOTH hard
argmax and soft probability assignments. Within supported cells, the pooled
all-galaxy count defines one target distribution shared by all environments
and both assignment methods. Each galaxy contributes its class weight times
the target cell fraction divided by that class's total weight in the cell.
Thus each environment integrates to one and has exactly the same joint cell
mix. This is coarsened distribution matching, not exact individual matching.
No red/blue preselection or outcome-dependent weighting is applied.

The primary estimand is the observed colour distribution of each inferred
class standardized to this common mass/redshift population. The fine check
uses 0.05 dex × 0.005 in z with its own supported population and pooled target;
it changes both resolution and support, so it is a sensitivity check rather
than a paired difference test. Original and common-support raw means are
both archived to distinguish support loss from reweighting.

## Results

Positive contrasts mean clusters are redder than voids. Intervals below are
conditional 95% sky-block bootstrap intervals for the mean colour contrast.

| Redshift | Original N | Primary retained N | Raw hard Δ g−r | Matched hard Δ g−r [95% interval] | Matched soft Δ g−r [95% interval] |
|---|---:|---:|---:|---|---|
| 0.15–0.20 | 534,070 | 533,037 (99.81%) | +0.04595 | +0.03724 [+0.03418, +0.04017] | +0.03110 [+0.02879, +0.03335] |
| 0.20–0.22 | 232,024 | 231,388 (99.73%) | +0.04645 | +0.03740 [+0.03328, +0.04106] | +0.03155 [+0.02859, +0.03422] |
| 0.20–0.30 | 611,197 | 607,579 (99.41%) | +0.04312 | +0.03742 [+0.03495, +0.03950] | +0.03114 [+0.02930, +0.03272] |

Matching reduces the hard-class mean contrast from 0.043–0.046 mag to about
0.037 mag; probability weighting gives about 0.031 mag. Distributions remain
strongly overlapping, but the colour contrast persists under both assignment
methods and both matching grids. Finer hard-class contrasts change by less
than 0.0002 mag. This supports a small descriptive conditional association,
not cleanly separated populations, new bimodality, or an environmental cause.
The two higher-redshift selections overlap and are not independent replications.

## Balance, uncertainty, and verification

All 48 environment histograms were independently rebuilt with direct weighted
histograms and reconciled to saved results. Achieved joint-cell target errors
are below 3e-14. Across primary hard-class means, the maximum pairwise mass
imbalance is 0.00125 dex and redshift imbalance 0.000122; on the finer grid these
fall below 0.000357 dex and 0.000047. These mean diagnostics supplement exact
cell balance; they do not prove identical continuous distributions within cells.
Primary hard-class effective sample sizes remain at least 12,331.

All 256 paired HEALPix nside=8 sky-block bootstrap draws are valid for every
window, matching grid and assignment method. Weights are refitted within each
draw to the fixed original target on fixed supported cells. Bands for each
bin's difference from the pooled matched distribution are pointwise 95%, not
simultaneous bands or a global distribution-equality test. Mean-contrast
intervals exclude zero in all twelve specifications. These are conditional
on fixed selection, mass/redshift estimates and environment predictions;
SED errors, completeness uncertainty, survey systematics and model-transfer
calibration are not propagated. No causal claim follows from matching alone.

All finite colour values enter normalization through explicit underflow and
overflow bins. Plotting zooms to 0.3–1.9 mag; no colour cut was applied. Histograms
use 1/30-mag bins. Source counts and unique TARGETIDs match all original windows.

Prepared analysis source: GraphWeb_DESI commit 7518c14. CPU allocation 59873074
on nid200063 completed analysis and independent verification with exit zero.
A rendering-only revision adds shared limits covering every peak, label space
and symmetric residual axes; original figures remain on Scratch. V2 numeric
CSV is byte-identical to the original. Primary hard and finer soft figures
were visually inspected. Small artifacts were copied with SHA-256 verification.

Scripts: `workflows/presentation/loa_matched_colour_20261016/`.
Evidence and CSV: `docs/evidence/loa_matched_colour_20261016/`.
Full output: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_matched_colour_20261016_v1/`.
`MANIFEST.json` records original output paths; `PRESENTATION_V2.json` records
final figure hashes. Archived figures use the final V2 content.
