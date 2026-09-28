# Published luminosity-function diagnostic — 2026-09-28

Executed `workflows/p12a_vac/compare_published_lf.py`; inputs and numerical
results are pinned in `evidence/p12a_published_lf_20260928/`. No catalogue,
selection, model weights or production recipe changed.

## Independent reference and comparison

Downloaded all four published total r-band North/South Vmax/SWML tables from
https://astro.dur.ac.uk/~cole/BGS/LFs.html (Moore et al. 2026;
https://arxiv.org/html/2511.01803). The recovered generator files contain log10
**cumulative** densities: verified in the archived luminosity_function.py.
Converted them into exact 0.25-mag bin averages before comparison to published
**differential** densities. Reference labels are z=.1, M-5logh and
h^3 Mpc^-3 mag^-1. This aligns nominal units, not all photometric/K/E conventions.

Generator wsys / published LF, ranges over both hemispheres and estimators:

| Magnitude interval | Ratio |
|---|---:|
| -21 to -20 | 0.918–0.958 |
| -22 to -21 | 0.875–0.937 |
| -23 to -22 | 0.794–1.015 |
| -24 to -23 | 0.097–0.217 |

Both wsys/nowsys were evaluated (JSON). Executed historical branch remains
unresolved. The recovered input LF has a much steeper extreme bright tail than
the published reference. This is not yet a measured LF of the generated parent;
it cannot by itself explain why v1 raw counts exceed completeness-corrected Loa.
Photometric definitions, realised HOD population and evolution must be separated.
No formal significance: covariance/systematics and matched photometry are absent.

## Controlled evolution response

Published global r Q=.78 versus recovered v1 Q=.67 changes apparent magnitude
by only -.044 at z=.5, holding all else fixed. Yet at a hypothetical reference
Mlimit=-23 the wsys cumulative density increases 1.509 times. This is a
**threshold response**, not a prediction of the survey n(z): actual limiting M
varies with colour/K-correction and redshift. Colour-specific published Q values
(.23 red,1.59 blue) imply opposite shifts (+.176,-.368mag at z=.5).
Their total-LF response values in JSON illustrate steepness only; they are not
red/blue number-count predictions. Changing Q without rederiving LF and colour
assignment is not a consistent replacement model.

## Interpretation and remaining discriminating work

1. Measure the realised v1 parent LF in volume-complete luminosity/redshift cells
   using stored M, then reconstructed M under both original and published K/E.
   Compare with both generator branches; retain central/satellite and colour bins.
2. Replay Q and K separately on the same parent rows, with fainter parent support
   explicitly checked, and compare joint z/r/colour distributions to Loa.
3. Derive independent Loa luminosities under published calibration; avoid calling
   the existing common-v1-K reconstruction an independent LF measurement.
4. Keep observational forward replay separate. Batch59013807 remains pending
   Resources at this check; no assignment result exists from that job yet.

Validation: monotonic input magnitude and cumulative-density grids; published
0.25-mag spacing; all interpolation inside table domain; SHA256 of every input.
No repair, density matching, retraining or release qualification was performed.
