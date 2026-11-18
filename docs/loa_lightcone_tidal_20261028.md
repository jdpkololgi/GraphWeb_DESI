# Loa full-cap tidal wedge display

Adds threshold-controlled T-web layers to the existing owner-private light-cone
viewer, project `appgprj_6ae18b46013c819183a767bbd43887f2`, without changing
the density/galaxy data, high-declination geometry, inference or training.

## Scientific contract

Uses the existing expanded U-Net seed17 four shared-state density draws, each
already targeting Gaussian R7 at native 6.766 Mpc/h spacing. For each cap
separately, apply the cuboid-integrated free-space potential Hessian over the
entire saved 3D field extent. The source is `(rho-1)*support` in prediction-covered
cells with positive fractional random support. All unmodelled exterior,
including the opposite cap, has delta=0. There are no periodic images, added
smoothing, fitted exterior or box-mean subtraction. This is a disclosed
closure-dependent exploratory display, not a calibrated full-survey tidal release.
Fractional footprint occupancy is not an angular completeness correction.

Classify the three 3D eigenvalues by strict `lambda > lambda_th` before
projection. Classes 0/1/2/3 mean void/wall/filament/knot. Thresholds -0.20 through
1.00 in 0.05 steps are encoded exactly: each uint8 stores the number of thresholds
strictly below one eigenvalue; 255 is missing. This is not eigenvalue rounding:
every encoded class is checked against its float64 calculation at all 25
thresholds. Other thresholds are not supported by this display encoding.
Float32 sampled eigenvalues are retained separately in Scratch HDF5 files.

Angular samples and masks match the published density grid exactly. Combine
class indicators with solid-angle weights through declination; optionally
average those fractions across four draws, never average eigenvalues or numeric
class labels before classification. Majority requires >50% (with 1e-10 numerical
tolerance); ties and no-majority pixels are mixed. Separate class-fraction layers
show composition. Fractions are projected occupancy/draw fractions, not calibrated
probabilities. Any incomplete angular support stays blank. Projection overlap
and unknown exterior tides preclude a 3D connectivity or recovery claim.

## Source and evidence

- Exporter: `workflows/presentation/loa_lightcone_tidal_20261028.py`.
- Operator: Illustris `workflows/sbi/e2e_footprint_tidal_operator.py` (unchanged).
- Output: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_lightcone_tidal_20261028_v1`.
- Viewer source: `/global/u2/d/dkololgi/loa_lightcone_site_20261028`.
- `manifest.json` binds source field, density manifest, exporter and operator
  hashes plus compressed asset hashes. Original field hashes match the existing
  density manifest. No source field was written.

Direct-sum error 4.44e-16; changed-padding error 6.66e-16; cube self-term error
1.11e-16. Full-cap trace errors <=2.09e-14 (NGC), <=1.43e-14 (SGC).
All eight draws completed; threshold classification and support agreement pass.
Compressed tidal assets total 28,828,705 bytes. Static site is about 253 MB,
below the 256 MiB deployment limit, without thinning or replacing old data.

Synthetic JS tests cover strict threshold equality, all four classes,
classification-before-average, angular weighting, mixed pixels and missing
support. Existing cone-geometry tests pass unchanged. Real-data browser QA
covers both caps, low/high declinations, thresholds -0.2/0.2/1.0, reference
pixels, fraction normalization, density preservation, unchanged galaxy counts,
hover, PNG export, polar geometry and 390px layout. No browser errors or external
requests. Class-map screenshots were inspected. This is local browser evidence;
native deployment status separately establishes publication, not scientific validity.

First browser check exposed a fraction of 1.0000000000000002; the renderer now
clips only numerical roundoff after rejecting deviations exceeding 1e-10. A
subsequent hover test targeted a canvas outside the viewport after clicking
export; the test now scrolls it into view and selects a supported in-wedge pixel.
Final browser QA passes. These were not scientific-data corrections.

## Execution

User approved one CPU node, 64 CPUs, desi/interactive, Scratch, at most 60 minutes,
with release when finished. Allocation 60428515, nid200012. Export step 0
completed in about 159 seconds after setup, with no new inference. Browser
steps 1/2/3 retain the failures described above; step 4 passed. Publication,
graph refresh and allocation cleanup are recorded below when confirmed.

Allocation released after 11m05s; `sacct` reports allocation COMPLETED 0:0,
export step 0 COMPLETED 0:0 (3m33s including setup), final browser step 4
COMPLETED 0:0. AST graph refresh step 5 failed before extraction with the
pre-existing `fcntl.flock` OSError 524. No bypass or global refresh from stale
output was attempted. Source synchronization required a credential renewal;
no access policy or credentials were written to source files.

## Publication

Owner-private deployment succeeded on 2026-10-28 at
https://loa-lightcone-slices.astrodak99.chatgpt.site . Exact source commit
`09910ac08ee0f4c34ce36faa1523fe9216f41545`; deployment
`appgdep_6ae1a1cc3e748191aa7ea41cd3d35c51`; saved version
`appgprj_6ae18b46013c819183a767bbd43887f2~appgver_d64c79163b98819183bf0ffd8e061d36`.
The matching archive is `site-deployment.tar` under the output directory above.
Local browser checks cover that source; a separate authenticated production
browser inspection was not run. The native deployment returned `succeeded`.
