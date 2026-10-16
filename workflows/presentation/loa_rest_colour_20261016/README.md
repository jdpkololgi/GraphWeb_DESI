# Observed/rest-frame colour and environment audit

Prepared for the user's 2026-10-16 request. Run only in an explicitly approved
CPU allocation; outputs use a fresh Scratch directory. `run.py` prepares the
join and selection, then computes trends with six worker processes. `render.py`
exports diagrams and trends after the analysis. Large FITS remain on Scratch.

Rest colour is **SDSS 0.1(g-r)**, `ABSMAG01_SDSS_G - ABSMAG01_SDSS_R` from Loa
FastSpecFit v1.0 main/bright. The reference redshift is 0.1, not zero. This is
not internally de-reddened. CIGALE masses/SFRs and P12-A environment probabilities
are retained. The input CIGALE VAC itself has no rest-frame luminosities.

The numeric cuts 0.75 < colour < 1.5 and colour < 0.75 are retained to expose
classification changes; 0.70/0.80 are predeclared sensitivity boundaries.
None is assumed to be a validated quenched/star-forming division. Quenched
fraction still uses CIGALE sSFR<1e-11/yr. Rest colours and CIGALE properties
share underlying observations and are not independent measurements.

1. Match the full Bright parent at 0.15<=z<0.30 by TARGETID. Verify unique IDs,
   FITS extension alignment, position <=1 arcsec and both FastSpecFit/Redrock
   redshift agreement <=0.001 normalized. Require finite magnitudes, positive
   g/r inverse variance, positive finite continuum chi-square and good redshift.
   Zero inverse variance flags synthetic-only magnitudes, excluded here.
2. Recompute faintest-20%/95th-percentile mass limits for each rest-colour split,
   retaining the original all/quenched/star-forming floor. Use the maximum across
   all three sensitivity boundaries and redshift bins, rounded up to 0.01 dex.
   This is empirical flux completeness conditional on available valid colours.
3. Preserve original observed diagrams; show hard argmax and probability-weighted
   environment counts. New rest diagrams use audited rest cuts. Paired observed/
   rest comparisons use identical valid-match rows and max(old,rest) mass limit.
   Audit missing/invalid matches by environment, mass and redshift. Missingness
   is not corrected by the flux-completeness estimator.
4. Retain 0.1/0.05 dex mass, 0.01 redshift, Kish>=20 controls within each colour;
   mass panel controls redshift only. Both soft-only and joint soft/hard support
   are recorded; unsupported comparisons stay unavailable. Use 256 sky-block
   bootstrap draws, conditional on fixed selection/colour/environment estimates.
5. Verify direct means, order-statistic percentiles, fresh source rows, FITS
   roundtrips and independent histogram bins. Save hashes, counts and manifests.

Interpretation limits: observed Legacy Survey and rest SDSS bands differ, so
changes measure frame + filter-system + fitting effects, not a pure K-correction
experiment. A redshift trend in these flux-selected samples can also reflect
mass threshold, internal attenuation, population evolution and selection.
No environment recalibration, inference retraining or causal claim is made.

Official source documentation:
- https://fastspecfit.readthedocs.io/en/latest/loa.html
- https://fastspecfit.readthedocs.io/en/3.0.0/fastspec.html
- https://fastspecfit.readthedocs.io/en/2.5.2/algorithms.html

Loa uses FastSpecFit 3.1.4; local v1.0 FITS headers were checked for the exact
SPECPHOT column names. These references define the product and magnitude
conventions; local counts and science results must come from the run records.

A separate blue-fraction diagnostic holds one common mass floor across every
0.01-wide redshift bin and compares both colour definitions on identical rows.
Its denominator includes objects outside the red/blue cuts; it is descriptive,
not a volume-corrected evolutionary measurement.

Execution completed on approved CPU allocation 59858561 (nid200021). Syntax
checks, four contract tests, 384 fresh FastSpecFit source rows, all selected
FITS roundtrips and independent original-VAC checks passed. The numerical
audit checked 1,364 reported curve vectors; every available calculation had
256 valid bootstrap draws. 46 unsupported controls remain explicitly unavailable.
Final plots are in `presentation_v2`: this expands sensitivity figures to all
12 predeclared combinations; the first render is retained on Scratch.
See `docs/loa_rest_colour_20261016.md` for results and release verification.
