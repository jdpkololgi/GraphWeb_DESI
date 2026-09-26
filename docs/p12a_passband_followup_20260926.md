# P12-A assignment and same-galaxy photometry follow-up — 2026-09-26

Diagnostic only: no catalogue repair, magnitude shift, refit, retraining or VAC
replacement. Follow-up to p12a_photometry_recipe_20260926.md.

## What a production record means

The earlier phrase “production receipt” meant a reproducibility record for the
specific catalogue: executed command, code revision, configuration/tables, input
versions and execution log. It is not a special DESI product or bureaucratic
requirement. Here it is needed to establish the origin of the measured
-0.8(z-0.1) magnitude relation and the effective assignment/mask implementation.

## Module-swap explanation now reproduced

A module swap asks Lmod to replace the currently loaded software module. The
archived scripts do not stop if this fails. In a separate shell, using the
retained 20240425-2.2.0 module directory identified in mock output headers:

- Loading fiberassign/main succeeded.
- Swapping to fiberassign/4.0.0 failed with status1: unknown module.
- Swapping to fiberassign/5.0.0 likewise failed with status1.
- After either failure, command lookup still returned the same
  `20240425-2.2.0/code/fiberassign/main/bin/fba_run`.

No assignment was executed. Only the child shell's module environment changed.
This directly reproduces the proposed failure mechanism. It does not prove the
historical failure: versions may have been removed later, or execution may have
used a local override/rerun. No historical stderr was recovered.

The retained Git history adds a useful bound: version string5.7.2.dev3588 was
introduced on2024-02-13 and changed on2025-04-29. It is compatible with the main
software used during that period, but is not a unique commit identifier. Thus
five matching version headers do not prove identical source revisions.
The five output files have retained modification times on2024-10-29–31, within
that version-string interval. This is consistent timing evidence, not an
execution timestamp guaranteed against later copies or timestamp preservation.
See module_swap_probe.txt and fiberassign_version_history.txt under the evidence
directory. Reproduction command: bash -lc 'bash workflows/p12a_vac/probe_archived_module_swap.sh'.
A module failure can affect assignment, but cannot by itself explain a deficit
already present in the raw bright parent.

## Output-flux convention caught and corrected in the audit

FastSpecFit3.1.5 io.py first divides input photometry by MW transmission;
fastspecfit.py lines99–105 then writes that parsed, already de-reddened flux
into its output METADATA table. Its output still includes MW_TRANSMISSION.
Therefore FastSpecFit output FLUX must not be divided by that column again.
This differs from the original Loa LSS catalogue convention.

The initial pilot made this double-correction in its selection and observed
imaging comparison. The final catalogue-parity check caught it. All affected
selection, imaging and SED-fit-quality diagnostics were rerun with the correct
convention; the tables here and final JSONs supersede the initial logs. Paired
synthetic filter offsets do not use observed flux or MW transmission, although
the corrected selection changes some sampled low-z objects. The original
production VAC did not use these FastSpecFit columns and is unaffected.
The corrected512-object check finds exact equality in both g/r between
FastSpecFit output flux and original Loa flux divided once by FastSpecFit's
transmission. Using each catalogue's own transmission changes g-r by a maximum
0.000197mag (median0.000083mag). All source hashes, unique matches and model-
colour replay checks pass; see CATALOGUE_PARITY.json and VALIDATION.json.
Logs same_galaxy_passbands_20260926.log, sdss_imaging_pairs_20260926.log and
passband_catalogue_parity_20260926.log contain the superseded first audit;
passband_corrected_convention_20260926.log is the corrected chain.

## Same-galaxy synthetic photometry

Read Loa FastSpecFit v1.0's main-bright hp00 catalogue, whose header records
FastSpecFit3.1.5 and speclite0.20. Selected512 existing VAC TARGETIDs in the
common southern footprint:32 randomly selected objects per dz=.025 bin from
0.15 to0.55, seed20260926, numerical12<=r<19.5. All512 have usable coefficients.
This is a stratified one-HEALPix pilot, not a survey-representative sample.

Reconstructed each fitted SED with the retained3.1.5 template implementation,
Chabrier2.0.0 templates, fitted coefficients, dust attenuation and velocity
dispersion. Applied SDSS2010 and DECam2014 g/r filters to that *same* SED.
Optical calculations omit infrared dust emission and IGM attenuation; saved
model-colour replay provides a numerical check on this approximation.

| z | Median SDSS r − DECam r | Median SDSS(g−r) − DECam(g−r) | N |
|---|---:|---:|---:|
| .15–.25 | +.093 | +.037 |128|
| .25–.35 | +.106 | +.055 |128|
| .35–.45 | +.132 | −.018 |128|
| .45–.55 | +.246 | −.089 |128|

Reconstructed observed DECam and rest-frame SDSS colours agree with the saved
FastSpecFit model colours to <0.00002mag. This verifies implementation parity;
it does not independently validate the fitted SEDs. Many objects have large
photometric reduced chi-square and model/imaging colour residuals. In the
high-z shell, restricting to22 objects with |observed−model g−r|<.03mag gives
r offset+.249 and colour offset−.089; the9 objects with RCHI2_PHOT<10 give
+.241 and−.070. These subset results support the sign/scale but are not unbiased
population estimates. See SED_ROBUSTNESS.json.

**Implication:** equal numerical SDSS/DECam r cuts need not select equal galaxies,
and the discrepancy can grow with redshift. At high z, DECam colour being redder
than SDSS also has the same direction as part of our DESI/mock colour offset.
However, these pilot medians cannot simply be subtracted from the full-sample
0.157mag offset: samples, model dependence and intrinsic populations differ.

Critically, the raw mock's measured -0.8(z-0.1) term may already embody a
convention or adjustment. Until its origin is known, interpreting R_MAG_APP as
an unadjusted SDSS magnitude or adding the offsets above would be unjustified.
The synthetic calculation does not reproduce the spatial Petrosian aperture.

For the same fitted SEDs, the median high-shell difference between synthesized
SDSS r K-correction and the mock GAMA table evaluated at the fitted rest colour
is only+.019mag (colour-K difference+.041mag). This pilot does not support
attributing the entire count deficit to a large r K-correction error alone.

## Direct SDSS imaging cross-check

Located the local DR16 mirror of SDSS dr13_final galaxy sweeps. Used the16 files
with the most nearby pilot targets among219 candidates, matched within1arcsec,
requiring a unique detection with SURVEY_PRIMARY (bit256), positive g/r fluxes;
no duplicate matches remained. The source definition is documented by
[SDSS resolve](https://www.sdss4.org/dr12/algorithms/resolve/), and flux units by
[SDSS magnitudes](https://www.sdss3.org/dr8/algorithms/magnitudes.php).

Found83 matches among512 pilot targets. That fraction is not completeness:
only16 files were read. Requiring Petrosian SNR>=5 in both g and r leaves60.
Extinction uses each catalogue's stored correction; differences in dust
convention have not been harmonized. No full clean-photometry flag selection
was applied; flags and per-object measurements are saved for inspection.

| z | N with SNR>=5 | SDSS Petrosian r − Legacy r | SDSS model r − Legacy r | SDSS Petrosian r − model r |
|---|---:|---:|---:|---:|
| .15–.25 |21|+.129|+.099|+.044|
| .25–.35 |20|+.265|+.204|+.047|
| .35–.45 |13|+.138|+.164|+.012|
| .45–.55 |6|+.435|+.302|+.076|

Each entry is a separate median, so columns need not subtract exactly.
These observations independently show that numerical magnitudes cannot be
assumed interchangeable. They are noisy, spatially restricted and selection-
conditioned, especially the six high-z objects; they do not establish a smooth
redshift correction or isolate passband, aperture, sky subtraction and dust.

## Remaining diagnostic work

1. Obtain the raw-catalogue producer's definition and executed transformation
   for R_MAG_ABS/R_MAG_APP, including the -0.8 term and archive changes. This
   determines whether the apparent mock magnitude already includes a mapping.
2. Extend the direct imaging comparison to representative common-sky targets,
   with consistent extinction, clean photometry, flux uncertainties and explicit
   surface-brightness/size dependence. Separately measure filter versus aperture
   effects; never infer them from equal column names or fit an arbitrary count offset.
3. Use the established physical magnitude definitions to replay parent selection
   and compare intrinsic luminosity/colour distributions and retention versus z.
4. Recover original FA stdout/stderr/configuration if available. The failed-swap
   explanation is now reproducible but historical attribution remains conditional.

Evidence: docs/evidence/p12a_passband_followup_20260926/. Scripts:
same_galaxy_passbands.py, sdss_imaging_pairs.py, passband_catalogue_parity.py,
probe_archived_module_swap.sh under workflows/p12a_vac. Compute allocation58894134, node nid001000, released after successful validation;
logs in sbi-logs. Existing E2E allocation58892711 was not modified or used.
