# Loa parent audit and empirical mass-complete subsets

This workflow reads the Loa CIGALE catalogue and the registered BGS Bright LSS
parent, then selects rows from the existing P12-A VAC. It does not retrain a
model, repeat environment inference, modify source catalogues or submit jobs.
Payload stages require an allocated compute node. Large outputs go to
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_mass_complete_20261016_v1/`.

## Selection

The parent audit uses positive finite CG15 masses, GALAXY spectra and
DELTACHI2 > 40 at 0 < z < 0.55. All-targeting and BGS Bright distributions are
separate. The Bright parent additionally requires ZWARN = 0 and extinction-
corrected r < 19.5 (South) or 19.54 (North). Repeated IDs in the eligible
parents are excluded rather than arbitrarily resolved. ID matches must also
agree within 1 arcsec and |delta z|/(1+z) <= 0.001.

For each 0.01-wide redshift bin, rescale each galaxy's mass to its flux limit:
`logM_lim = logM + 0.4 (r - r_limit)`. Take the 95th percentile for the faintest
20% separately in the combined, quiescent and star-forming populations, and
adopt the largest of the three. Each population must contain at least 100
galaxies; unsupported bins remain unavailable. This is an empirical flux-
completeness estimate conditional on the observed, fitted sample, not a
guarantee of absolute survey completeness.

Two sub-VACs preserve every original P12-A column and add selection/provenance
columns. `zdependent` applies the bin-specific limit over 0.15 <= z < 0.55.
`fixed020030` applies one threshold over 0.20 <= z < 0.30, equal to the maximum
bin limit in that interval rounded up to 0.01 dex. Both require SUPPORTED and
an unambiguous CIGALE galaxy match, but do not require valid SFR, age, dust or
AGN measurements. A fresh exact ID/sky/redshift/mass check binds the joins.

`QUENCHED_SSFR11` is 1 for SFR/M < 1e-11 per year, 0 otherwise, and -1 when
SFR is invalid. Zero SFR is valid and quenched; negative/nonfinite SFR is not.
All masses use the source CG15 solar-mass convention without an additional
h conversion. Added `SOURCE_ROW` is the zero-based original VAC row.

## Relationships

Primary results use the fixed 0.20-0.30 sample. Alternatives include the
redshift-dependent selection, the flux-limited parent, an additional 0.1 dex
mass margin, and the full 0.15-0.55 redshift-dependent sub-VAC. Property-specific
validity masks avoid imposing age/dust/AGN availability on the quenching
sample. Soft P_WEB-weighted and hard argmax classes share a reference mass/z
distribution on cells with Kish effective count >= 20 in every class/method.
Mass bins are 0.1 dex, checked against 0.05 dex; redshift bins are 0.01. The
stellar-mass panel controls redshift only. Results report retained common
support and redshift coverage; they do not impute unsupported cells.

Quenching thresholds -10.5, -11 and -11.5 in log10(sSFR/yr^-1) are compared.
The primary choice is -11. Observed g-r remains a descriptive property; it
is not treated as a rest-frame UVJ classifier. Mean log SFR/sSFR use a labelled
-14 floor for zeros; the quenched flag uses the unfloored ratio.

Intervals are 95% from 256 paired HEALPix nside=8 sky-block resamples, with
fixed control cells and reference weights. They exclude mass/SFR errors,
uncertainty in the estimated completeness boundary, model transfer and
correlated environment-posterior uncertainty. Different selection samples
can have different retained support, so differences between their curves
are not an isolated causal effect of the mass cut. Survey selection, fibre
assignment and surface brightness can still bias a mass-selected sample.

## Execution and validation

Run `run.py --stage prepare`, then `--stage analyse`, `--stage analyse_full`,
and `--stage render` in the approved allocation using cosmic_env Python with
PYTHONPATH/PYTHONHOME/PYTHONUSERBASE/LD_PRELOAD unset and PYTHONNOUSERSITE=1.
Stages refuse to overwrite JSON/FITS results. Use a new versioned output
directory for a new scientific selection; do not delete an old release.

`test_run.py` contains small selection and weighting arithmetic tests.
`validate.py` independently reconstructs limiting masses and checks every
original column against fresh source reads at 257 rows in each output.
`validate_relationships.py` independently reconstructs weighted quenching
curves, raw means on identical rows and contrasts for all three thresholds
and both mass-bin widths. Run the catalogue validators on compute nodes.

Method references:

- [Empirical limiting-mass construction](https://www.aanda.org/articles/aa/pdf/2010/16/aa12801-09.pdf).
- [Quenching definitions and measurement systematics](https://academic.oup.com/mnras/article/485/4/4817/5376508).
- [DESI DR2 luminosity-function sample and redshift quality](https://lss.fnal.gov/archive/2025/pub/fermilab-pub-25-0845-ppd.pdf).
