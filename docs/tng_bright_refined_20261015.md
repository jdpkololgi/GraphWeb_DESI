# Refined z=0 TNG photometry and BGS selection proxies

The remaining quenched-fraction excess survives actual DECam filter integration,
BGS colour/fiber cuts under two profile proxies, and finer mass bins. Its size
is larger than in the preceding single-colour approximation. Dust and genuine
stellar spectra remain unavailable, so this is not a validated survey mock.

## Why mass selection did not already account for brightness

Matching stellar mass and observing redshift matches the population distribution
in those two variables. Brightness also depends on luminosity per unit stellar
mass, which varies with stellar population, star formation, metallicity and dust.
At fixed mass and redshift, a flux limit can preferentially exclude quenched
objects. A high mass threshold can reduce this bias if it lies safely above the
completeness limit for all relevant populations. The earlier logM≥10.5 cut does
not establish that condition over z=.2-.3 (the previous empirical completeness
curve rises to roughly10.9). Neither that threshold nor mass-bin reweighting
makes the TNG luminosity selection equivalent to DESI's.

## Results

Each row below uses identical Loa cell weights and common support across its
selection arms. Offsets average the four equally weighted environment quarters.
Different rows have slightly different common support; the bin-width check is
therefore a sensitivity test, not a perfectly fixed-population contrast.

| Mass bins | Population | No brightness cut | Previous colour approximation | DECam r only | BGS n=1 proxy | BGS n=4 proxy |
|---|---|---:|---:|---:|---:|---:|
| 0.25 dex | Shared support | 29.15 pp | 4.62 pp | 8.62 pp | 8.63 pp | 8.62 pp |
| 0.25 dex | 10.5≤logM<11.5 | 28.99 pp | 7.23 pp | 12.18 pp | 12.20 pp | 12.18 pp |
| 0.10 dex | Shared support | 30.09 pp | 4.51 pp | 8.69 pp | 8.70 pp | 8.69 pp |
| 0.10 dex | 10.5≤logM<11.5 | 29.60 pp | 6.97 pp | 12.04 pp | 12.06 pp | 12.04 pp |

Nominal n=1 residual95% spatial intervals are7.98–9.33 and11.37–13.05pp
for the coarse panels;8.05–9.31 and11.22–12.83pp for the fine panels.
They exclude photometric/SFR/mass modelling uncertainty and transfer error.
At high mass/coarse bins, selected TNG quarters69.2,70.2,70.7,75.6% compare
with Loa55.3,58.1,59.9,63.6%. Their environmental contrasts6.44 versus8.30pp
are also not identical. Fine-bin contrasts5.10 versus7.86pp.

![Refined selection overlay](/global/u2/d/dkololgi/fmsu_tokyo_20261010/tng_bright_refined_20261015/refined_overlay.png)

The actual DECam r response is redder than the previous nominal SDSS-r
wavelength approximation. The combined continuum/filter update makes predicted
r magnitudes about0.083mag brighter at the median in this sample, allowing
additional galaxies into the selection. The inferred quenched excess increases.
This is a combined change, not a separately identified filter-only effect.

At z=.205,14,517/66,691 TNG galaxies pass integrated-r selection and14,510
pass the n=1 proxy; at z=.295 the corresponding counts are3,162 and3,151.
The n=4 proxy removes none beyond the flux cut here. Therefore the tested fiber
cuts do not explain the residual. This does not establish that real imaging,
blending, surface-brightness or spectroscopic selection is negligible.

## What improved and what remains approximate

1. Fresh source extraction of dust-free SDSS griz magnitudes and stellar
   3D half-mass radii from the existing snapshot99 group catalogues, aligned
   to the same66,691 held-out Subfind IDs. No new cosmological data downloaded.
2. Positive piecewise power-law f_nu continuum with four fitted flux nodes.
   Iterative adjustment reconstructs the input SDSS2010-integrated griz
   magnitudes; then redshift this continuum and integrate actual DECam2014
   g/r/z response curves from the installed speclite1.0.0 data. Planck18
   luminosity distances; observing bin midpoints .205,...,.295. This replaces
   the monochromatic g-r K approximation. It is **not** a fitted stellar
   population spectrum: griz cannot resolve spectral breaks, emission lines
   or dust. Observed g additionally relies on extrapolation below rest g.
3. Apply installed `desitarget.cuts.isBGS_colors` southern Bright logic:
   broad g-r/r-z colour bounds,12≤r<19.5, the fiber-magnitude criterion,
   and fiber-total bright-star safeguard. Independently replay the applicable
   inequalities to check exact selection parity. No BGS Faint colour cut.
4. Fiber flux proxy integrates a circular Sérsic n=1 or n=4 profile through
   a1.5arcsec diameter aperture convolved with1arcsec FWHM Gaussian seeing.
   Set projected half-light radius to .75 times the catalogue 3D stellar
   half-mass radius. The conversion and profiles are assumptions, not
   measured light profiles. Fiber-total flux equals isolated model fiber
   flux; blends, contaminants and star/galaxy separation are not simulated.
5. Same Loa CG15 properties and explicit Bright eligibility as before.
   Recompute rank quarters within each selection and mass–redshift cell;
   use TNG R7 truth and Loa P12-A posterior-mean R7 trace. No inference or
   graph recomputation. Matching is binned and does not match SFR by design.
6. Repeat0.25dex and0.10dex mass bins with dz=.01, ≥80 rows per cell in
   every arm. Coarse full/high support88,971/68,918 Loa rows; fine87,581/
   68,849 out of89,948/68,918 eligible. Mass capped at11.5 in all panels.
7. Same200 spatial resamples across TNG observing-redshift placements;
   independent Loa sky-block resampling. Keep support and reference weights
   fixed. One rare replicate loses a small fine cell; discard invalid
   replicates for the affected interval, leaving at least199/200 per curve.

No internal attenuation model, stellar-particle spectra, true epoch evolution,
photometric noise, image quality/masks, Gaia separation, fiber assignment or
spectroscopic success were added. Realistic dust cannot be recovered uniquely
from these four dust-free bands. The failed public-package download was replaced
with read-only use of existing DESI filters; no project installation was changed.

## Conclusion

Brightness selection still substantially reduces the excess beyond mass matching.
A residual of roughly9pp overall/12pp high-mass is stable to these mass-bin and
profile-proxy checks. The earlier4.7/7.2pp values were not stable to improved
filter treatment. Neither set establishes robustness to realistic dust and
stellar spectra, or assigns the residual to cosmic evolution/SFR measurement.
The next decisive improvement requires a dust-attenuated catalogue or spectra
for these same z=0 galaxies; a new redshift snapshot is not required for that
selection/measurement step. Actual evolution remains a separate later test.

## Validation and artifacts

Source: `workflows/presentation/tng_bright_refined_20261015.py`; independent
validator: `workflows/presentation/validate_tng_bright_refined_20261015.py`.
Outputs: `~/fmsu_tokyo_20261010/tng_bright_refined_20261015/`: PNG/PDF overlay,
`RESULTS.json`, photometry inputs/predictions and validation receipts.

1260 direct cell-mean checks;3 prior high-mass curve replays;288 independent
speclite f_lambda integrations using5A instead of10A sampling (maximum
magnitude difference6.04e-7); source photometry/size parity; direct selection
inequality parity; monotonic fiber fractions and normalized reference weights.
Input-band reconstruction max error1.82e-10mag is numerical fit precision,
not evidence that the reconstructed spectrum is physically correct.
Figure inspected. Source/filter/cut hashes retained. Graphify unavailable.

References: [TNG field specifications](https://www.tng-project.org/data/docs/specifications/),
[speclite filters](https://speclite.readthedocs.io/en/v0.15/filters.html),
[DESI filter mapping](https://desispec.readthedocs.io/en/latest/_modules/desispec/io/filters.html),
[DESI BGS selection paper](https://arxiv.org/abs/2208.08512).
