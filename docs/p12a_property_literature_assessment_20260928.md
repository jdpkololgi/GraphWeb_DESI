# Interpreting modest galaxy-property contrasts

Literature assessment, 2026-09-28. Overlapping mass–sSFR distributions are
compatible with published cosmic-web results, but do not establish that our
inferences are accurate. No matched quantitative comparison with another finder
or the old model has yet been performed. No new compute, fitting, inference,
mock replacement or repair was performed for this review.

## Primary sources

Sources accessed 2026-09-28; different observables, samples and environment
definitions prevent a pooled numerical effect-size comparison.

| Study | Method and relevant result |
|---|---|
| [Alpaslan et al. 2015, GAMA](https://arxiv.org/abs/1505.05518) | Mass-matched large-scale environment samples have similar colours, energy outputs, luminosities and morphologies. Removing mass control amplifies differences. Direct precedent for overlapping distributions, not a CIGALE sSFR benchmark. |
| [Eardley et al. 2015, GAMA](https://arxiv.org/abs/1412.2141) | Tidal-tensor void/sheet/filament/knot luminosity-function differences are explained by local density. Uses 4 and 10 h^-1 Mpc smoothing; relevant to our scale, but luminosity is not sSFR. |
| [Alam et al. 2019, SDSS](https://academic.oup.com/mnras/article/483/4/4501/5255205) | Overdensity delta8 and tidal anisotropy alpha5: quenching depends mainly on density, with main trends reproduced by a halo-quenching model. Additional direct large-scale tidal effects are weak. This does not imply raw cluster–field contrasts vanish. |
| [O’Kane et al. 2024, SDSS](https://arxiv.org/html/2409.09028v1) | Projected DisPerSE: filament–field star-formation differences persist after mass matching but disappear after also matching local density. Their ~1 Mpc filament-distance definition differs from ours. |
| [Hoosain et al. 2024, RESOLVE/ECO](https://arxiv.org/abs/2401.09114) | DisPerSE colour/gas gradients survive mass control but are dominated by group and halo composition; an additional small gas effect remains in low-mass haloes. |
| [Kraljic et al. 2018, GAMA](https://academic.oup.com/mnras/article/474/1/547/4430643) | 3D DisPerSE detects passive-fraction and star-forming colour/sSFR gradients at fixed mass. Density-test interpretation depends on scale; Appendix C shows near-erasure when matching small-scale DTFE density. |
| [Darvish et al. 2017, COSMOS](https://www.ipac.caltech.edu/publication/2017ApJ...837...16D) | Stronger low-z SFR changes occur, particularly for satellites; star-forming-only trends differ from total-population trends. These medians are not our mean log-sSFR estimand. |
| [Nandi & Pandey 2025, preprint version read](https://arxiv.org/html/2507.18614v1) | SDSS DR18 tidal classification, 8 Mpc smoothing and zero threshold: finds mass-matched quenched-fraction differences. Counterexample to a blanket claim of no residual web dependence; thresholds and selection differ from ours. |

[Libeskind et al. 2018](https://arxiv.org/abs/1705.03021) compares twelve finders
and documents method dependence alongside shared features. Our R=7 h^-1 Mpc,
threshold0.2 tidal knot is not equivalent to a bound cluster or a narrow
DisPerSE filament. Mpc and Mpc/h smoothing scales must also be distinguished.

## What our plots actually show

The same galaxies enter every panel with different P_WEB weights; panels are
normalized by class probability weight and displayed with a logarithmic colour
scale. They are not disjoint populations. Both modes can remain visible while
their relative weights differ. This is an interpretation of plotting mechanics,
not a fitted mixture decomposition or proof of identical mode locations.

Existing checksum-bound [wedge results](p12a_cigale_environment_20260928.md):

| Comparison | Knot–void mean log sSFR | Knot–void low-sSFR fraction |
|---|---:|---:|
| Raw CG15 probability weighted | -0.469 dex | +16.40 percentage points |
| Mass/redshift controlled, probability weighted | -0.205 dex | +7.77 points |
| Controlled argmax | -0.270 dex | +10.23 points |
| Controlled argmax, max P>=0.8 | -0.344 dex | +14.06 points |
| Controlled CG5 probability weighted | -0.186 dex | +7.93 points |

Supports differ: 94,028 primary, 91,902 argmax, 37,260 confident and94,027 CG5.
Stronger hard/confident contrasts are consistent with membership dilution but
do not isolate it, because the samples change. The adjusted low-sSFR fraction
changes from43.48% to51.25%; it is not flat. Spatial-bootstrap intervals exclude
zero but omit property and transfer systematics.

We control stellar mass and redshift, not local density, halo mass or satellite
status. Thus this is not a pure geometric-tide effect. Even calibrated P_WEB
does not generally make sum(P_k*y)/sum(P_k) an unbiased true-class property mean.
Calibration, class discrimination and preservation of property contrasts are
different questions. Stronger astrophysical contrast is not a training objective
or a valid criterion for choosing the best web finder.

## What remains unresolved

The [halo48 mock review](../../TNG/Illustris/docs/p12a_halo48_readiness_20260925.md)
reports useful coverage and knot Brier0.03420 versus0.03423 on matched draws.
That small relative change is not a prevalence-baseline comparison. The
[conditional/golden audit](../../TNG/Illustris/docs/p12a_golden_conditional_results_20260925.md)
explicitly limits its qualification scope. Marginal coverage and observer parity
do not alone establish accurate DESI class separation or contrast retention.

Proposed discriminating checks, not executed in this literature review:

1. Use already-exposed existing mock validation products to compare true classes
   with saved predictions at the same smoothing, threshold and selection. Report
   reliability, confusion and per-class proper scores against prevalence and
   density baselines. Do not open reserved phases or replace mock populations.
2. Test attenuation of known environment–property relations. If mocks lack
   validated SFR/mass truth, inject explicitly synthetic properties with null
   and known-amplitude relations. This tests sensitivity, not galaxy physics.
3. Compare an independently constructed density/filament reference on identical
   wedge galaxies and properties. Match mass/redshift support and distinguish
   scale/definition changes. Agreement is corroboration, not matter truth.
4. Repeat soft/hard/confident comparisons on common support, tighten control
   bins, use control-bin permutations, and compare old/new models on identical
   rows. Add density and group controls separately; they answer an incremental
   geometry question. Standardized sSFR CDF/difference plots would expose small
   contrasts more directly than broad 2D logarithmic heatmaps.

Weak truth contrasts recovered faithfully would favour a physical explanation;
strong truth contrasts suppressed by predictions would implicate inference
attenuation. Good mock recovery but poor matched DESI agreement would flag
transfer/selection. Literature alone resolves none of these alternatives.
Keep the existing provisional model and all real-DESI coverage/mismatch caveats;
the instruction to stop mock repair remains in force.
