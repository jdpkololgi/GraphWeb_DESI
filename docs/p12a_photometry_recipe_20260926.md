# P12-A photometry and generator-provenance investigation — 2026-09-26

Diagnostic only. No repair, changed selection, refit, retraining or VAC replacement.
The existing VAC remains provisional. Only already exposed ph006 was read.

## Findings

### A numerical identification of the colour-conversion recipe

Recovered the author's public Abacus generator and tables:
https://github.com/amjsmith/shared_code/tree/a40a5d3be570a6ce9d5a0698268b04864d4d3af2/abacus_mocks

A deterministic raw-row sample (every 320th row; 200,184 rows, of which 65,402
have 0.15<=z<0.55 and r<19.5) reproduces
`G_R_OBS = G_R_REST + K_g(Z,G_R_REST) - K_r(Z,G_R_REST)`.
The largest colour residual in these science rows is 1.84e-7 mag. This identifies
both tables, cubic interpolation in rest colour, clipping to coefficient-table
colour limits, and use of observed/RSD Z rather than Z_COSMO. Below z=0.5,
linear rather than cubic colour interpolation leaves RMS residuals 0.003–0.006
mag; using Z_COSMO leaves measurable scatter. Above z=0.5, the code uses linear
redshift extrapolation based on the z=0.4–0.5 secant, with linear colour
interpolation of the secant coefficients.

This is strong numerical evidence for that mapping, not a recovered production
commit or proof of the luminosity/colour population-generation recipe. The
public code describes SDSS Petrosian photometry and GAMA K-corrections; equality
to DESI's Legacy Surveys flux-based selection has not been demonstrated.

### A previously unpinned r-band magnitude convention

With a native-like flat LCDM distance in Mpc/h, the stored columns obey

`R_MAG_APP - R_MAG_ABS - DM(z) - K_r(z,c) ~= -0.8*(z-0.1)`.

Fitted slope is -0.800183, intercept -1.99e-5 mag; after fitting the line,
RMS residual is 5.79e-7 mag. Fixing the coefficient to -0.8 leaves RMS 5.02e-5
mag (maximum 1.03e-4); the distance calculation neglects radiation and is not
an exact replay of the generator cosmology. The r residual before accounting
for this term is about -0.08 mag in the first shell and -0.33 mag at z>0.5.
The colour identity above does not depend on a distance model.

The public cut-sky function directly converts its evolved absolute magnitude;
its LF parameter table has P=1.8 and Q=0.7. It does not by itself document the
extra relation seen in the stored columns. The measured -0.8 term must not be
identified with the LF's Q=0.7 parameter. It could encode a stored absolute-
magnitude evolution convention or a subsequent apparent-magnitude adjustment;
we have not established which column was transformed or its motivation. It is
not evidence to apply an offset to DESI or to the mocks. A naive absolute-
magnitude comparison that ignores this relation would be wrong.

### Archived revisions are demonstrably different

Independent 10,000-row samples of ph006 v0, v0.1/old and current v0.1 all
reproduce the same colour conversion to float precision. Their r-band mappings
differ: the fitted term is approximately -0.8(z-0.1) for v0, -1.6(z-0.1)
for v0.1/old and -0.8(z-0.1) for current v0.1. Using Z_COSMO in the distance
reduces the older-revision residual scatter substantially, whereas current
v0.1 matches using Z. Older residuals remain 0.0015–0.0019 mag even with
Z_COSMO, so this is not an exact old-distance replay. The additional old
magnitude term and changed distance mapping are real version-provenance leads;
calling them a double correction or assigning a motivation would be premature.
These are separate row samples, not a joined-galaxy or comparable count test.
The current P12 ph006 source is current v0.1, not the old directory. No archived
revision was substituted. See ARCHIVE_FINGERPRINT.json and its reproducible
photometry_archive_fingerprint.py script.

### The colour discrepancy survives matching apparent magnitude and redshift

Compared the actual trained-on ph006 observed catalogue against Loa on the
existing common angular support, SGC only and Loa PHOTSYS=S. Applied identical
*numerical* 12<=r<19.5 limits. Within cells of dz=0.01 and dr=0.1, require at
least 20 objects from each catalogue, then standardize mock cell colour
histograms to DESI cell counts for this diagnostic. No catalogue is altered.
This controls z/r composition and removes mixed northern/southern real imaging;
it does not assume that the physical magnitude definitions are equivalent.

| Redshift | DESI median g-r | Standardized mock median | DESI minus mock | DESI retained in common cells |
|---|---:|---:|---:|---:|
| 0.15–0.25 | 0.990 | 1.005 | -0.016 | 99.86% |
| 0.25–0.35 | 1.301 | 1.328 | -0.027 | 99.71% |
| 0.35–0.45 | 1.603 | 1.584 | +0.020 | 99.04% |
| 0.45–0.55 | 1.736 | 1.579 | +0.157 | 92.88% |

Colour histograms use 0.02-mag bins; these are descriptive histogram medians,
not confidence intervals. High-shell common cells retain 38,304/41,242 DESI
and 21,002/21,246 mock objects. The remaining discrepancy therefore cannot be
explained solely by different coarse-shell z/r mixtures. Spectroscopic
selection and physically different passbands/magnitude estimators remain
confounded with colour-population modelling.

The z>0.5 branch is not an established cause of the deficit. Relative to
continuing the polynomial (not an alternative truth model), its median r
change is -0.021 mag for 167 sampled bright objects at z=0.50–0.55: it makes r
brighter in that comparison. The mismatch already develops below z=0.5.
Changing this branch without an empirical SED/passband comparison is unjustified.

### Fibre-assignment version discrepancy is systematic in sampled receipts

Five date-spaced first-tile examples cover 2021-05-14, 2022-10-20, 2023-03-25,
2023-10-09 and 2024-03-23. Scripts request 4.0.0 (first) or 5.0.0 (others);
all output headers record **5.7.2.dev3588**. Thus this is not one isolated file.

The scripts source the main DESI environment, then request a module swap and
run fba_run, without fail-fast handling. Current installed LSS's wrapper
captures subprocess output without checking its return code at that call.
A failed swap could therefore leave an already available executable in use.
This is a viable mechanism, not proof of the historical event: no stdout/stderr
receipts were found in the bounded altmtl6 log/output-file inventory. Stale
scripts, local overrides and reruns are also possible. The header records the
loaded implementation; the archived shell request does not certify that version.
Current source is explanatory evidence, not a historical production receipt.

This FA uncertainty cannot alone explain the upstream count deficit: the raw
bright parent already has fewer high-z objects than successful DESI on the same
sky, before assignment and redshift-quality losses.

## Remaining cause-identification work

1. Recover the raw v0.1 producer's receipt defining R_MAG_ABS, the measured
   -0.8(z-0.1) relation, exact LF/colour tables and any post-generation changes.
   Numerical fingerprints narrow the request considerably but do not replace it.
2. Compare DESI synthetic SDSS and Legacy Surveys photometry for the *same*
   high-z targets, with an explicit Petrosian-versus-model magnitude treatment.
   This separates passband/flux-estimator effects from the mock colour/LF model.
3. With that mapping established, compare the evolving absolute-luminosity and
   rest-colour distributions and stage retention. Do not fit an arbitrary
   offset to the observed count ratio.
4. Recover executed FA/catalogue software and mask receipts. Preserve actual
   output-header versions; do not claim module requests are executed versions.

## Reproducibility

Evidence: `docs/evidence/p12a_photometry_recipe_20260926/` with source tables,
SHA256 hashes, FINGERPRINT.json, RAW_SAMPLE.npz, MATCHED_CELLS.json/npz and
ASSIGNMENT_RECEIPTS.json. Source scripts are in workflows/p12a_vac:
photometry_recipe_fingerprint.py, photometry_matched_cells.py and
assignment_receipt_audit.py. Compute allocation 58892691, node nid001016;
logs in sbi-logs/photometry_*_20260926.log. No external downloaded code was
executed by the numerical fingerprint; the table evaluation is independent.
