# Matched observed-colour distributions (completed)

Retain original observed-colour/CIGALE-mass samples and frozen P_WEB. Match the
joint mass/redshift distribution across hard and soft inferred classes on
shared support: minimum effective sample size 20 in every class under BOTH
methods. Target is the pooled galaxy distribution on supported cells, held
fixed across methods and bootstrap draws. No colour selection is applied.
Primary cells 0.10 dex × 0.01 z; finer sensitivity 0.05 dex × 0.005 z.
Matching is exact at cell resolution, not at every continuous mass/redshift.
Report achieved cell weights, within-cell residual mean mass/z differences,
effective sample sizes, exclusions and original versus matched colour means.

The 256 paired HEALPix nside=8 sky-block resamples refit per-cell weights to
the same target; support is held fixed. At least 205 valid draws required.
Residuals use the equal-environment pooled matched distribution; their bands
are pointwise, not a global equality test. Mean cluster-minus-void intervals
are descriptive and conditional on selection, fixed masses/redshifts and
fixed environment predictions; they exclude inference/SED systematic errors.
All finite colours enter normalization including explicit outlier bins.

Three window workers, eight BLAS threads each. Intended budget: one CPU node,
64 CPUs, desi / interactive, 30 minutes, Scratch, then release. Approved allocation 59873074 completed on nid200063; analysis and verification passed. Inputs are the frozen observed_original FITS products from
`loa_rest_colour_20261016_v1`; new Scratch destination
`loa_matched_colour_20261016_v1` is created exclusively.
Syntax and two synthetic contracts pass: known mass-mixture differences are
removed, outlier probability retained, and unsupported cells excluded.
