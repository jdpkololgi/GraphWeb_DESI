# CIGALE properties, independent density, and eigenvalue geometry

Computed 2026-10-14. The independent point-count estimator reproduces the directions and approximate amplitudes of all seven current CIGALE trends. The individual-galaxy correlations are modest for both independent density and inferred trace. This supports the environmental trend measurement on this selected sample; it does not determine a universal expected amplitude.

## Direct matched comparison

Current 15-band CIGALE values are joined by exact TARGETID to the previously verified 94,415-row neighbour-count cache. All property panels require the same finite, physically valid values. There are 89,809 eligible interior galaxies; 84,275 remain in common mass/redshift support. Stellar mass uses the 89,809-row redshift-only comparison. Other properties use identical 0.1-dex stellar-mass / 0.01-redshift strata and identical reference weights for every estimator. The minimum effective sample size is 20 per group per stratum.

Independent density is the self-excluded neighbour count inside a 7 h^-1 Mpc top-hat sphere (10.34585 Mpc), divided by the frozen expected tracer count. The original calculation used 316,480 surrounding full-VAC tracers, including points outside the displayed wedge. Galaxies must pass the stored boundary-distance cut. Raw-count quartiles provide a second non-neural comparison without the expected-density normalization. All quartiles are defined within dz=0.01, keeping tied values together.

| Property | Density Q4-Q1 | Raw-count Q4-Q1 | Trace Q4-Q1 | Soft cluster-void |
|---|---:|---:|---:|---:|
| Mean log stellar mass [dex; z controlled] | +0.1287 | +0.1305 | +0.1055 | +0.1055 |
| Mean log SFR [dex] | -0.2067 | -0.2114 | -0.1713 | -0.1903 |
| Mean log sSFR [dex] | -0.2077 | -0.2125 | -0.1721 | -0.1908 |
| Quenched fraction [percentage points] | +7.6904 | +7.8881 | +6.4059 | +7.3643 |
| Mean log fitted AGE [dex] | +0.0103 | +0.0109 | +0.0086 | +0.0098 |
| Mean AV [mag] | -0.0529 | -0.0529 | -0.0413 | -0.0486 |
| Mean FAGN [fraction] | +0.0064 | +0.0062 | +0.0040 | +0.0056 |

Density produces a 7.69-percentage-point quenched-fraction rise, with a 16-84% sky-block interval of 7.29-8.03 points. Soft web classes give 7.36 points (6.91-7.67); inferred-trace quartiles give 6.41 (6.10-6.73). The paired density-minus-soft difference is 0.33 points (0.05-0.68). The curves are close, although not identical estimands: quartiles and probability-weighted classes select different environmental mixtures. These are endpoint contrasts, not correlation coefficients.

The mean log-sSFR density contrast is -0.208 dex, corresponding to a ratio of geometric means of about 0.62. On the same common rows, ordinary Spearman correlations of log sSFR with density and inferred trace are -0.146 and -0.120. Centering within-cell ranks gives -0.081 and -0.069. Thus modest individual-galaxy association coexists with a clear ensemble shift. The narrow uncertainty on a mean does not imply a tight individual-galaxy relation.

All 128 paired sky-block bootstrap replicates were valid. Intervals describe spatial sampling with fixed strata, quartile thresholds and reference weights; they omit correlated SED-fit uncertainty and environment-transfer uncertainty. AGE is the supplied SED-fit parameter; it has not been identified as a mass-weighted stellar age. FAGN is a fitted contribution fraction, not an AGN incidence or a spectroscopic activity classification. Its small positive shift should not be compared directly with literature AGN-host fractions.

The earlier 9.29-versus-9.18-point density/soft comparison used the older property cache and different common support (65,566 galaxies). The new 7.69-versus-7.36 result uses current CIGALE values and the present seven-property support; the frozen environment model has not changed.

## What the literature supports

- Kauffmann et al. (2004), https://arxiv.org/abs/astro-ph/0402030: in a complete SDSS sample, density measured in 2-Mpc projected cylinders with a +/-500 km/s window has a strong sSFR association at fixed mass. Below 3e10 solar masses the low-to-high-density median sSFR change exceeds a factor of ten, and the strongest dependence is on sub-Mpc scales. This is evidence against declaring every weak environment trend expected. Their median, density extremes, tracer selection, mass range and aperture differ from ours.
- Balogh et al. (2004), https://arxiv.org/abs/astro-ph/0311379: projected fifth-neighbour density and Gaussian densities on 1.1/5.5-Mpc scales show changing star-forming fractions, while the H-alpha equivalent-width distribution within the star-forming population changes little. A changing mixture can therefore coexist with weaker within-population trends.
- Peng et al. (2010), https://arxiv.org/abs/1003.4747: red fractions exhibit approximately separable mass and environmental dependences. This motivates matched mass/environment comparisons, but their colour-based populations do not prescribe a numerical value for our CIGALE sSFR threshold.
- Alam et al. (2019), https://arxiv.org/abs/1801.04878: colour-based quenching is primarily associated with spherical density delta8, with subtler anisotropy alpha5 dependence at fixed density; their halo model accounts well for the main trends. Weak additional large-scale tidal dependence must not be confused with weak total quenching-density dependence.
- Kraljic et al. (2018), https://arxiv.org/abs/1710.02676: GAMA shows mass/type gradients with distance to filaments and nodes and lower star-forming-galaxy sSFR nearer filaments at fixed mass. Their residual geometry result makes geometry beyond density worth measuring; it does not imply every tidal summary should show a large marginal trend.

Assessment: the signs and density-dominated ordering are consistent with these studies. A 0.2-dex average sSFR shift and roughly eight-point quenched-fraction change are plausible for this broad, selected 7-h^-1-Mpc comparison. The literature does not establish that exact amplitude as expected for this DESI/CIGALE sample. The strongest present evidence is the same-galaxy independent-density agreement. A quantitative literature benchmark would require matched tracer selection, aperture, stellar-mass support, quenching definition and summary statistic.

## Non-neural environment measurements to add

| Observable | Added information | Main requirement |
|---|---|---|
| Projected fifth-neighbour Sigma5 = 5/(pi d5^2), with a fixed velocity window | Adaptive local crowding; tests whether small-scale trends are stronger | Tracer completeness, edge correction, documented velocity window |
| Counts at multiple radii; inner/shell density ratio | Compact group-like overdensities versus broad overdensities | Same tracer population; density and redshift controls |
| Neighbour-direction quadrupole or shape tensor | Angular anisotropy at fixed neighbour count | Poisson correction, mask tests, redshift-distortion sensitivity |
| Distance to a non-neural filament skeleton or nodes | Geometrical proximity beyond density | Reconstruction persistence/smoothing and edge tests; not a direct tidal tensor |
| Group richness, projected group-centric radius, close-pair separation | Halo-scale crowding and interaction-scale environment | Group-catalogue membership uncertainties |
| Density-grid Poisson reconstruction and traceless tidal shear | Closest physics-based comparison to inferred tides | Selection/mask treatment, smoothing, galaxy bias and boundary assumptions |

The first and multiscale-count tests follow established approaches in Balogh et al. above; filament distances follow Kraljic et al. A classical tidal reconstruction is described by Alonso et al. (2016), https://academic.oup.com/mnras/article/460/1/256/2608903. The direction-quadrupole suggestion is a geometric diagnostic, not a claim that the local point-cloud shape equals the gravitational tidal tensor. Compare its property associations at fixed density, mass and redshift. None of these additional estimators was computed in this bounded cached analysis.

## Why the eigenvalue contours are asymmetric and diagonal

For this catalogue lambda1 <= lambda2 <= lambda3. A plot of a larger eigenvalue on the vertical axis is restricted to y >= x. Reflection symmetry about y=x is therefore not expected. Exact equality of two eigenvalues represents a degeneracy; only equality of all three makes the tensor isotropic.

Write T = (t/3) I + S, with t = tr(T) and tr(S)=0, so lambda_i = t/3 + s_i. Changing t at fixed S moves every eigenvalue equally and hence moves a pair parallel to y=x. The difference lambda_j-lambda_i removes this common shift and measures one component of anisotropy. Positive covariance and variation in both trace and shape produce tilted contours; the tail is also affected by evolved density statistics and galaxy selection. Ordering alone creates the one-sided boundary, but does not uniquely explain the full contour shape.

The existing full-range histogram bins give DESI pair correlations of about 0.693, 0.632 and 0.845 for (1,2), (1,3) and (2,3). Abacus truth gives 0.685, 0.617 and 0.863. The diagonal structure is already present in truth. These approximate binned moments describe the displayed population mixtures, not the posterior covariance of one galaxy.

The companion figure re-expresses the (lambda1,lambda3) histogram as pair midpoint versus gap. The midpoint is not pure density: it mixes trace and traceless shape; a true density/shape decomposition uses all three eigenvalues. The display is rebinned from saved histograms, with no new inference or posterior sampling.

With the usual normalized Poisson convention, tr(T)=delta_m for the identically smoothed matter field. Point counts instead estimate a selected galaxy overdensity delta_g; the large-scale relation delta_g approximately b delta_m is an approximation. Our top-hat tracer count is therefore a useful independent proxy, not an exact measurement of the model trace.

## Reproduction and validation

Scripts: workflows/presentation/desi_density_properties_20261014.py and desi_corner_geometry_20261014.py. User explicitly requested login-node execution; all numerical libraries were restricted to one thread. No Slurm job submitted or modified. The analysis reuses cached counts rather than rebuilding a spatial index.

Verified density-cache SHA256, source receipt row count, exact extracted TARGETIDs, current P_WEB/redshift parity and eigenvalue ordering. The 5.4-million-row property catalogue was not fully rehashed: its frozen receipt hash is inherited and that limitation is explicit in provenance.json. Extracted rows have their own SHA256. Synthetic shared-reference, paired-identity, tie and within-cell-rank checks passed. Figures were visually inspected. Source and all deliverables are hashed in MANIFEST.json.

The wedge sample was inherited from an older positive-mass/SFR cache; replacing property values does not remove that historical selection. Counts share survey selection with the ML inputs even though no neural weights enter their construction. Results are not full-footprint or volume-complete estimates.
