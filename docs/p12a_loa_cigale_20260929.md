# CIGALE properties matched to the provisional Loa posterior VAC

The user-supplied DR2 CIGALE catalogue is joined by TARGETID to the frozen
P12-A full-survey VAC. Neither the Loa selection nor any posterior is changed.
The enriched product retains every original row and column and adds CIGALE
properties, metadata and match flags. It is a property-augmented provisional
VAC, not a newly calibrated environmental catalogue.

## Inputs and matching

Source: `/global/cfs/cdirs/desicollab/users/zouhu/vac/dr2/dr2_galaxy_sedfitting_v1.0.fits`.
The 33,346,076-row catalogue includes multiple observations of some TARGETIDs.
CG15 combines five Tractor bands and ten spectrophotometric bands; CG5 uses
five Tractor bands. The supplied catalogue description says AGN contributions
are modelled and cautions that low-redshift SFR can be affected by Tractor
modelling. Header units are used: stellar mass in solar masses, SFR in solar
masses/year, age in Myr and attenuation in magnitudes. No dedicated CIGALE
fit-quality flag was identified; the spectroscopic CHI2 column is not used as
a SED-fit quality measure.

Of 5,436,413 VAC rows, 5,436,352 have a TARGETID present in the source. After
requiring angular separation <1 arcsec and |dz|/(1+z_VAC)<0.001, 5,019,523 have
exactly one consistent CIGALE row. All accepted coordinates agree exactly.
There are 416,135 TARGETIDs with more than one consistent source row: these
remain in the VAC with CIGALE_AMBIGUOUS=true, CIGALE_MATCHED=false and missing
property values. This conservative exclusion is not evidence those galaxies
have poor fits. Resolving their observation/coadd provenance could recover
additional properties; selecting the first or highest-SNR fit is not assumed
scientifically equivalent to the Loa observation. Another 755 VAC rows have
no acceptable unique match and are not marked ambiguous. Matching rejected
2,357 source rows for coordinate/redshift inconsistency; this is a row count,
not an additional number of missing VAC galaxies.

`CIGALE_SOURCE_ROW` is zero-based in the source FITS table; -1 means no accepted
match. Separate CIGALE_Z/SEP_ARCSEC/DZ_NORM retain join diagnostics. The original
VAC Z is never replaced. Properties for unmatched rows are NaN, not zero.

## Analysis

Plots require an accepted unique match, SUPPORTED=true, CIGALE SPECTYPE=GALAXY,
and finite positive stellar mass and SFR. Nonpositive SFR is not silently
classified as quenched. No new AGN exclusion or unvalidated SED-quality cut is
made. The CG15/CG5 comparison requires both fits to be usable on the same rows.

Environment weights are the existing P_WEB in void, sheet/wall, filament, knot
order. Every galaxy contributes to all four probability-weighted panels;
these are descriptive relationships with inferred probabilities, not unbiased
estimates of true-class galaxy populations. Two-dimensional histograms share
axes, density normalization and colour limits.

Full-survey control cells use dz=0.025 and CG15 log-mass width 0.1 dex over
8<=logM<13; every class must have effective N>=20 in a cell. Mass is standardized
to a common redshift distribution. sSFR and the descriptive low-sSFR fraction
(log10(sSFR/yr^-1)<-11) use a common mass/redshift distribution. Both fit variants
use the CG15 mass cells and identical galaxies/reference mix, so their comparison
is not an independent control on each fit's own stellar mass. Intervals are
16–84% from 64 nside=8 sky-block bootstrap replicates, holding the reference mix
fixed; these exclude SED parameter covariance, model-transfer and other systematics.

The original wedge (RA 120–160 deg, Dec 14.5–30.6 deg, 0.20<=z<0.30) is also
reproduced with the established plotting workflow, 0.25-dex mass cells and 96
spatial bootstrap replicates. It uses this new catalogue and conservative
matching, not the old goodPhoto-filtered cache; exact sample counts need not
match the earlier figures. Its historical CG5 sensitivity panel uses CG5's own
mass and potentially different support, unlike the new controlled comparisons.

## Qualification

Real-DESI environmental coverage remains unverified. There is a known training
mock population mismatch above z~0.35, severe at 0.45<z<0.55. Redshift trends can
reflect selection/model mismatch and must not be interpreted as physical
evolution without independent validation. Lower-redshift results are not
certified automatically. Probability weighting can attenuate environmental
contrasts. The newly available SED catalogue does not resolve these limitations.

## Results

The analysis has 4,989,187 usable supported CG15 galaxies, including 2,327,331
at 0.20<=z<0.30 across the survey footprint and 102,634 in the original wedge.
“Full sky” in the figure labels means the full available Loa survey footprint,
not all-sky completeness. There are 166 QSO-classified unique matches; these
are retained in the enriched catalogue but excluded from the galaxy plots.

Controlled knot-minus-void differences:

| Sample / fit | Mean log mass (dex), z controlled | Mean log sSFR (dex), mass/z controlled | Low-sSFR fraction difference |
|---|---:|---:|---:|
| Original wedge, CG15 | +0.113 | -0.197 | +7.54 percentage points |
| Full footprint, 0.20–0.30, CG15 | +0.105 | -0.186 | +7.13 pp |
| Full footprint, 0.20–0.30, CG5 | +0.107 | -0.180 | +7.67 pp |
| Full selected range, CG15 | +0.104 | -0.188 | +7.34 pp |
| Full selected range, CG5 | +0.106 | -0.180 | +7.79 pp |

For the full-footprint 0.20–0.30 CG15 comparison, the 16–84% sky-bootstrap
interval is [-0.189,-0.184] dex in the sSFR contrast and [7.03,7.22] pp in the
low-sSFR fraction contrast. Low-sSFR fraction rises from 44.98% (void-weighted)
to 52.10% (knot-weighted) after standardization. CG5 changes the absolute sSFR
normalization by about 0.19 dex, but gives a similar environmental contrast.
These are modest, resolved descriptive differences amid strongly overlapping
bimodal distributions; their small sampling errors are not a measure of total
astrophysical or environmental-inference uncertainty.

The matching-completeness figure is essential: conservative duplicate exclusion
is strongly redshift dependent. Around z=0.50–0.55 only roughly 46–52% of the
supported class probability weight survives into usable CIGALE properties,
compared with about 97% at 0.20–0.30. This is principally a match-selection
restriction, not absence of TARGETIDs in the supplied catalogue. It is also
environment dependent. The full-range and high-z plots are therefore
selected-subset diagnostics and cannot establish a full-population redshift
trend, even apart from the known mock-transfer mismatch. Redshift/mass control
does not prove missing properties are ignorable.

## Products and reproduction

The full enriched catalogue is
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits`.
Large FITS and parquet products remain on scratch. The similarly named JOIN
FITS is an intermediate only; FULL_VAC is the deliverable with all original
posterior columns. JOIN.json and PRODUCT.json record input/output hashes,
row checks and exact original-column preservation.

[Full-survey figures (8 pages)](figures/p12a_loa_cigale_20260929/FULL_SURVEY_CIGALE_environment.pdf)
and [reproduced wedge figures (4 pages)](figures/p12a_loa_cigale_20260929/original_wedge/LOA_CIGALE_environment.pdf).
Numerical results, source header/units and SHA256 manifests are alongside them.
Final plot run: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_plots_20260929_v3`.

Run on a CPU allocation using cosmic_env, from GraphWeb_DESI:

```bash
python workflows/p12a_vac/loa_cigale_full.py --out NEW_JOIN_DIR
python workflows/p12a_vac/loa_cigale_product.py --root NEW_JOIN_DIR
python workflows/p12a_vac/loa_cigale_plots.py --join NEW_JOIN_DIR/DESI_LOA_P12A_CIGALE_FULL_VAC.fits --out NEW_PLOT_DIR
python workflows/p12a_vac/property_environment.py --test
```

CPU allocation 59073991 executed the scans and plotting. The first plot attempt
hit a FITS big-endian/Parquet conversion incompatibility; explicit native-endian
conversion fixes it. The final plot run also uses float64 probability reductions
and records the completeness curves numerically. No fitting or inference ran.
