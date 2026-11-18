# Loa cube-projection website

Six density/tidal/galaxy cube presets: NGC/SGC × three anchor redshift ranges
0.15–0.28, 0.28–0.41, 0.41–0.55. Expanded U-Net seed17, four existing shared-state
draws. These are presentation diagnostics, not new posterior samples or a
qualified DESI field release.

## Scientific contract

- Each cube is 32³ native cells at 6.766 Mpc/h, side 216.512 Mpc/h. UI permits
  16³/24³/32³ subcubes, three Cartesian projection axes and one-cell to full-cube
  depths. No added smoothing or survey-support display mask.
- Anchors are chosen deterministically near catalogue-coordinate medians,
  preferring support-rich, fully computed geometry. Density and environment
  are not used for selection. Anchor redshift bins do not further filter the
  galaxies inside the selected spatial cube.
- Density projection is the arithmetic mean of native density through the
  selected depth. Fixed logarithmic colour limits are display-only.
- Classify each 3D cell by the number of ordered tidal eigenvalues strictly
  greater than λ, then project equal-volume class fractions. The displayed
  most-frequent class is a depth summary, with ties explicitly grey; it is not
  a new 2D T-web classification. Four-draw mode averages cell/draw fractions,
  not a thresholded mean tensor or a claimed calibrated probability.
- Tidal arrays use the existing free-space cuboid operator, prediction-covered
  density sources of requested width 96 or 144 native cells, zero overdensity
  outside the finite source and in uncovered cells. No periodic wrapping.
  Source geometry and coverage fractions are recorded in QA.json. These two
  contexts probe sensitivity, not the true unknown cosmological exterior.
- Every parent-VAC galaxy in the exact half-open spatial cube/slab is overlaid;
  no row thinning, tile-edge, support or property-eligibility cuts. Coordinates
  come from the same validated parent catalogue used for field geometry. This
  is the selected 0.15–0.55 parent sample, not every raw Loa BGS row. Markers
  are not galaxy-environment classifications. Axes are equatorial Cartesian
  comoving Mpc/h, not angular sky or radial line-of-sight projections.

## Sources and verification

Generator: workflows/presentation/loa_cube_projections_20261028.py.
Scratch output: graphweb_desi/outputs/loa_cube_projections_20261028_v1.
Field and parent source hashes, operator/script hashes, cube bounds and marker
TARGETID hashes are bound in manifest.json. Native density is Float32; tidal
Float32 rounding errors are <1e-6. All finite/positive-density, eigenvalue-order,
trace, byte round-trip and fraction-sum checks pass. Threshold-near-boundary
classifications retain the disclosed Float32 numerical sensitivity.

Galaxies per cube (no thinning): NGC 62,956 / 28,334 / 4,589;
SGC 76,951 / 31,286 / 5,686. Context changes at λ=0.2 affect 0.45–1.99% of
cell/draw labels across these cubes, not a calibration or convergence result.

Independent JavaScript coordinate-ramp fixtures and Python reference pixels
check axes, cropping and projection arithmetic. Real-data checks cover each
cube with all axes, selected draws/four-draw means, threshold monotonicity,
fraction normalization and exact galaxy membership. Browser QA covers all
six cubes, both contexts, threshold changes, crop offsets, whole-cube and
one-cell projection, draw averaging, overlay visibility, fraction maps and
390px layout; no page errors or external requests.

Allocation: 60424442/nid200011, user-approved one node/32 CPUs/desi/interactive,
Scratch, up to 30 minutes; release after verification and packaging.
Site: appgprj_6abdae1cd14c8191bf07f24c6b0a6a59; preserve owner-only audience.
Page: cubes.html, linked from the homepage alongside tidal.html.
