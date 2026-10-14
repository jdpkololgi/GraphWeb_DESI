# Current-P12-A presentation figures and DESI movie (2026-10-13)

User requested a current-catalogue version of the supplied mass-binned quenched-fraction figure, a DESI environment animation, an updated eigenvalue corner plot and threshold sweep, then an explicit redshift-control diagnostic. One CPU interactive allocation (59745887, nid200140, 30-minute limit) was explicitly approved. No inference, posterior sampling, training or catalogue mutation is involved.

Sources: immutable P12-A halo48 20260924 v1, full Loa catalogue 20260925 v1, unique CIGALE join 20260929 v1; already-exposed ph006 mock evaluation only. All panels use 0.20 <= z < 0.30. The preparation scripts verify recorded source hashes and exact draw-derived class fractions. Real-DESI calibration and selection limitations remain unchanged.

## Mass and redshift

Four requested log10 stellar-mass bands: [9.8,10.4), [10.4,10.6), [10.6,11.0), [11.0,11.6). Quenched means log10(sSFR / yr^-1) < -11, using positive finite 15-band CIGALE mass/SFR, unique consistent galaxy matches and supported environments. P_WEB gives soft group weights. Within each broad mass band, standardize to a common 0.1-dex mass / 0.025-redshift cell distribution, requiring Kish effective count >=20 in every comparison group. The continuous axis bins sum(EIGENVALUE_MEAN); it is not a measured density truth field. Class and trace controls are fitted separately and may retain different supports (in this run, retained row counts agree in each mass band).

The follow-up uses identical galaxies for raw, common-redshift, and common-fine-mass/redshift estimates. Cluster-minus-void contrasts in percentage points:

| log mass | Raw | Redshift | Mass + redshift |
|---|---:|---:|---:|
|9.8–10.4|1.051|0.963|0.726|
|10.4–10.6|9.507|7.418|6.788|
|10.6–11.0|10.944|10.754|9.853|
|11.0–11.6|8.110|8.092|6.753|

Redshift balancing particularly attenuates the 10.4–10.6-band slope but does not explain the residual association in general. The lowest band is near the quenched-fraction floor; other controlled contrasts are 6.8–9.9 pp. The full 0–1 y-axis makes slopes appear shallower than the supplied screenshot's truncated axis. Mass dependence, soft-membership attenuation, flux-limited population selection, property errors and model transfer are possible contributors; the curves do not identify causal quenching or quantify transfer accuracy. Intervals use 96 nside8 sky-block bootstrap replicates at fixed reference mix, excluding property/transfer systematics.

## Movie

20 s, 1920x1080, 24 fps, H.264/yuv420p. Geometrically fixed RA[120,160), Dec[14.5,30.6), z[.20,.30) wedge; supported galaxies with no CIGALE or probability cut. Display 18,000 uniformly sampled galaxies without replacement, seed20261013. Planck18 comoving redshift-space Cartesian coordinates rotated into local east/north/line-of-sight. Colors are argmax P_WEB, opacity .20+.80 max(P_WEB). Camera rotation and introductory color fade only: no galaxy motion, cosmic-time evolution, posterior resampling or coherent-field claim. Task-local imageio-ffmpeg 0.6.0 supplies H.264 because the system encoder lacks libx264.

## Corner and threshold sweep

21,052 ph006 evaluation galaxies with saved natural_weight; 512 saved joint triplets each. DESI uses 128 uniformly sampled original spatial cores (seed20261013); 17 contribute supported z=.20–.30 galaxies, totaling 16,126 rows, 512 draws each. This is a bounded descriptive subset, not a full-survey threshold sweep. No P_WEB or eigenvalue-based selection. Saved draws exactly reproduce the selected catalogue probabilities at threshold .2; mock fractions reproduce the existing A07 comparison. Native observed and mock selection functions remain different.

Corner: dimensionless ordered eigenvalues. Main window [-.6,1.8] excludes at most 0.404% marginal probability; normalization includes these tails, and a full-range companion is supplied. 50/90% contours use weighted 2D histograms, one-bin Gaussian display smoothing, and total-probability levels. Gold denotes both plotted eigenvalues exceeding .2, not all three. Full-range companion spans all used draws. Threshold sweep uses exact ordered-eigenvalue empirical CDF differences, equivalent to counting eigenvalues strictly above each threshold in each saved joint triplet. Fractions sum to one. No multiplication of marginal probabilities.

## Artifacts and reproducibility

- Durable presentation package: `/global/homes/d/dkololgi/fmsu_tokyo_20261010/desi_visuals_20261013/`.
- Compute outputs: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/tokyo_visuals_20261013_v1/`.
- Entry points: `workflows/presentation/desi_mass_movie_20261013.py`, `desi_corner_sweep_20261013.py`, `desi_redshift_control_20261013.py`.
- Exact plotted arrays/statistics, source hashes and QA are packaged beside SVG/PDF/PNG/MP4 files.

Three estimator sanity checks pass: constant response, removal of composition-only differences, retention of a true within-cell contrast. Visual checks include 720p figure previews and four movie views; H.264 stream has 480 frames at 24 fps and 20 s duration. Final packaging receipt records hashes and scheduler completion. `graphify` is unavailable in this shell; targeted direct inspection was used and graph refresh could not run.

Final verification: allocation59745887 COMPLETED0:0 and released after19m21s. Final task-local ffmpeg full decode completed0:0. Initial system-ffmpeg encode/decode limitations were resolved with the task-local build. Four vector PDFs have embedded fonts and no raster image objects; all arrays, IDs, probabilities and packaged source hashes pass. Final numerical replay reproduces the reported values.
