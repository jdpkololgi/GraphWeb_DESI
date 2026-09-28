# DESI Loa provisional environmental VAC: visual atlas and release notes

Prepared 2026-09-28. This package uses the existing frozen U-PATCH/P12-A halo48
model, its existing training mocks and the original full-survey posterior values.
No new inference, retraining, mock replacement, selection change or posterior
correction was performed. The mock-repair investigation is stopped at the user's
request. Its running preparation job 59013807 was cancelled; existing outputs
and scientific evidence were retained.

## Status and mandatory interpretation

- **Provisional environmental inferences**, with real-DESI coverage unverified.
- A known training-mock population mismatch, especially above **z~0.35** and
  severe at **0.45<z<0.55**.
- Redshift-dependent posterior trends may reflect selection/model mismatch;
  they **must not be interpreted as physical evolution without independent
  validation**. The saved posteriors do not marginalize this mismatch.
- Lower-redshift results are **not automatically certified** by the smaller
  count discrepancy.

This is a documented provisional product, not a science-qualified release.
The atlas is a visualization of per-galaxy inferences, not a reconstruction of a
spatially joint posterior field. Morphological agreement is a useful descriptive
check, not independent truth or a coverage measurement. No external publication
or transfer is performed by this packaging step.

## Product and integrity

Original FITS (unchanged):
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_full_20260925_v1/DESI_LOA_P12A_HALO48_FULL_SURVEY_VAC.fits`

SHA256: `cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc`.

The full catalogue contains **5,436,413 unique TARGETIDs**; **5,404,568** are
supported and **31,845** retain explicitly flagged null posteriors. All rows
remain provisional. No rows were removed from the product. The current whole-file
hash matches the original completion receipt both before and after plotting.
Checks also pass for finite/ranged coordinates, probability bounds and sums,
ordered eigenvalues, monotonic quantiles, support flags and every historical
quality-bit count. This verifies technical integrity, not scientific calibration.

Machine-readable results, exact display geometry, column dtypes and counts:
[INTEGRITY_AND_VIEWS.json](figures/p12a_loa_atlas_20260928/INTEGRITY_AND_VIEWS.json).
The [original completion receipt](evidence/p12a_loa_full_20260925/FULL_VAC_COMPLETE.json)
and its neighbouring source manifest/plan preserve the production provenance.
The original `DRAW_INDEX.json` and posterior-draw shards remain alongside the
FITS on Scratch; this plotting run does not independently revalidate every draw
shard (the original completion records that validation).

## Atlas

[Download the three-page PDF](figures/p12a_loa_atlas_20260928/LOA_environment_atlas.pdf).

1. [Full-depth slice](figures/p12a_loa_atlas_20260928/01_cosmic_web_slice.png):
   48,798 galaxies in a fixed equatorial slab, RA 130–230 degrees, absolute
   equatorial height <10 comoving Mpc, 0.15<=z<0.55. The three panels show
   identical positions, the most probable environment and largest class
   probability. The 36 unsupported rows are grey crosses in inference panels.
   Curves mark z=0.35 and 0.45, not physically inferred structures.
2. [Close-up](figures/p12a_loa_atlas_20260928/02_cosmic_web_closeup.png):
   6,044 galaxies in a fixed transverse ±180 Mpc, radial-plane 650–950 Mpc
   window of the same slab (actual z=0.1523–0.2309). All are supported.
3. [Four probability maps](figures/p12a_loa_atlas_20260928/03_four_environment_probabilities.png):
   the same close-up shown separately for void, sheet, filament and knot.
   All four use identical 0–1 colour scales. Low probabilities remain visible.

Geometry was specified without inspecting inferred classes. No confidence cut,
class balancing, downsampling, spatial smoothing or connecting curves were used.
The display subset does not change catalogue selection. Coordinates are
Planck18 comoving **Mpc**, from observed redshifts, retaining redshift-space
distortions; this is not a real-space reconstruction. If r is the radial distance,
x=r cos(DEC) sin(RA−180°), y=r cos(DEC) cos(RA−180°), and slab height=r sin(DEC).
A 10,001-point distance table over z=0.15–0.55 is interpolated for plotting only.
The model sees three-dimensional context beyond this thin display slab, so
visual two-dimensional morphology need not correspond exactly to its classes.
Sky masks, radial selection and redshift-space distortions affect appearance.

## Column guide

| Columns | Meaning |
|---|---|
| TARGETID, RA, DEC, Z | Galaxy identity, ICRS angles in degrees, observed redshift |
| CAP, CORE_ID | Original production region and processing-core identifiers |
| SUPPORTED, QUALITY | Support boolean and bitwise quality mask; use the documented bit definitions |
| BASE_EIGENVALUES | Saved deterministic encoder eigenvalues; distinct from posterior summaries |
| EIGENVALUE_MEAN, Q05/Q16/Q50/Q84/Q95 (each prefixed EIGENVALUE_) | Three-component ordered eigenvalue posterior summaries |
| P_WEB | Four joint-draw-derived probabilities, ordered void, sheet, filament, knot |
| BOUNDARY_MPC | Original production boundary diagnostic in Mpc |
| NTILDE_MPC3 | Original fitted expected-density input, in Mpc^-3 |

Targets are ordered lambda1<=lambda2<=lambda3 at Gaussian R=7 Mpc/h and target
epoch z=0.2. At threshold lambda=0.2, class 0/1/2/3 counts eigenvalues above the
threshold within each joint draw. These are environments at galaxy positions:
void-class galaxies are not a map of empty volume. The atlas uses argmax(P_WEB)
only for display; probabilities and eigenvalue posteriors remain primary.
Maximum class probability describes model-conditional concentration, not
empirically established accuracy on DESI. Eigenvalues alone do not infer
filament orientations.

Quality bits (unchanged): 1 unsupported/null; 2 sparse z>=0.45; 4 boundary<R/h;
8 boundary<2R/h; 16 outside training feature envelope; 32 provisional;
64 unresolved selection mismatch in z>=0.35 shells. A zero mismatch bit below
z=0.35 is not certification. Do not silently redefine the released selection by
filtering quality bits; report any analysis subset explicitly.

## Reproduction and preservation

Canonical script: `workflows/p12a_vac/visual_atlas.py`. Run with cosmic_env on a
CPU allocation and an absent output directory:

```bash
python workflows/p12a_vac/visual_atlas.py \
  docs/evidence/p12a_loa_full_20260925/FULL_VAC_COMPLETE.json \
  docs/figures/p12a_loa_atlas_20260928
```

Validated on CPU allocation 59020269, nid004216; final layout replay on
CPU allocation 59020572. The report binds the script
hash, allocation and catalogue checksum; `ARTIFACT_SHA256.json` binds all atlas
outputs. The small atlas/evidence are committed here. The large original FITS
and draw shards remain on Scratch and require durable collaboration archival
before long-term distribution; this step has not copied them to shared storage.

Further work in this stream is provisional delivery/documentation, not mock
repair. Renewed training or scientific qualification would require a separate
explicit decision and validation programme.

## Galaxy-property follow-up

[Current VAC versus saved CIGALE wedge properties](p12a_cigale_environment_20260928.md)
provides raw and mass/redshift-controlled star-formation comparisons. This
property cross-check does not change the VAC or its provisional status.
