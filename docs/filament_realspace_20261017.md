# Real matter structure, finite sampling and redshift distortions in the same three regions

The original true-wall string contains resolved anisotropic matter structure in
real space. It is not created solely by sparse galaxy sampling or redshift
positions. However, the previously sheet-dominated example is strongly sensitive
to which galaxies enter the aperture and to redshift distortions. These are
three illustrative development regions, not a population estimate.

## Fixed regions and comparisons

Keep TARGETIDs6032112,3743465,1493453 from the previous examples. Use the same
R=7/0.6766=10.346 Mpc observer-coordinate radius, centered on each galaxy in its
respective coordinate system. Use Z_COSMO for real-space positions and observed
Z for redshift-space positions. The original observed-coordinate counts and
covariances reproduce exactly. Compare:

1. Native dense matter from the existing ph006 2048^3 TSC counts grid, built from
   all field+halo A+B subsamples (10% of particles), cell size0.9765625 Mpc/h.
   This is genuinely dense relative to galaxies, but not the full particle set
   or sub-cell structure. Use both total mass and positive excess over cosmic
   mean for covariance. No new matter simulation or tidal operator is fitted.
2. All available galaxies in the original public HOD CutSky, before the bright
   parent and observation selection. These are denser (46/29/25 inside R), not
   an infinitely dense or volume-complete galaxy census.
3. 512 uniform subsets without replacement, of17/14/13 galaxies respectively,
   matching each region's actual DESI-mock observed neighbour count. Every
   real/redshift pair uses identical galaxy IDs. Counts exclude the centre.
4. The actual observed galaxies, carrying exactly the original redshift-space
   R-sphere membership, at real versus redshift-space coordinates. Also retain
   a separate fixed-real-sphere membership control to diagnose aperture turnover.

Uniform thinning matches local observed counts, not the complete DESI selection
function or global nbar(z). It conditions on observed counts and does not include
count uncertainty. Magnitude, halo occupation, targeting and fibre selection
remain distinct from random thinning. Bright-parent thinning is also retained
where possible; regions2/3 have only13/12 bright parents inside the real sphere,
less than their14/13 observed redshift-sphere neighbours, so that specific control
is unavailable rather than sampling with replacement or silently changing N.

## Results

L is covariance linearity, not filament probability. L/P/S sum to one for every
individual shape. Medians of separate components need not sum to one.

| Region | Dense real matter L | All available real galaxies: N, L | Thinned real galaxies: median L | Same observed IDs: real L -> redshift L |
|---|---:|---|---:|---|
| True-wall string6032112 | .547 | 46, .687 | .686 | .664 -> .664 (N17) |
| Previously sheet-dominated true wall3743465 | .443 | 29, .765 | .741 | .293 -> .042 (N14) |
| True-filament error1493453 | .587 | 25, .565 | .565 | .599 -> .762 (N13) |

**Original true-wall string:** total matter L/P/S=.547/.216/.236, positive-excess
matter .646/.231/.123. The real matter and denser galaxies already contain
anisotropic structure. 503/512 (98.2%) real-space thinned galaxy draws are
line-dominated; L's central95% thinning range is[.454,.858]. Keeping exactly the
17 original observed IDs gives L=.6639 in real space and.6637 in redshift space;
the original string does not depend on redshift distortion of those IDs.
RSD can nevertheless amplify a different fixed cohort: the23 observed galaxies
inside the real-space sphere change L=.544 -> .728. This distinction matters:
no claim that velocities have zero effect on this whole region is justified.

**Previously sheet-dominated true wall:** total matter L/P/S=.443/.265/.291;
positive-excess matter .693/.157/.150. The denser real galaxy sample is strongly
line-dominated (L=.765), and504/512 (98.4%) thinned draws remain line-dominated,
L95% range[.432,.915]. The actual14 redshift-selected neighbours have real
L/P/S=.293/.416/.292, becoming.042/.644/.314 after RSD: P-L increases+.479.
Membership is also crucial. The real sphere has11 observed neighbours, versus14
in redshift space, with10 common,4 entering and1 leaving. Holding the11 real-sphere
IDs fixed yields L=.637 -> .704, still line-dominated. Thus the earlier “sheet”
geometry cannot be read as a straightforward image of a dense matter sheet:
observation selection, aperture turnover and coordinate distortion all matter.
This does not contradict its true tidal-wall label.

**True-filament prediction error:** dense matter L/P/S=.587/.053/.360 and
positive-excess matter .769/.020/.211. The denser real galaxies are already
line-dominated.422/512 (82.4%) thinned real draws are line-dominated, L95% range
[.202,.785]. The same13 observed IDs change L=.599 -> .762 under RSD; P-L changes
-.262. The fixed-real12-ID cohort likewise changes L=.565 -> .743. This case
contains real anisotropy that velocities enhance; P12-A's original wall argmax
remains a genuine error against the frozen tidal truth.

## What finite sampling can and cannot establish

An additional unbiased matter-count control samples512 sets from the interpolated
matter density, with replacement on the display quadrature grid. Line-dominance
frequencies are422/512,357/512,400/512 (82.4%,69.7%,78.1%). These are matter-sampling
controls, not a galaxy-bias prescription or DESI selection emulator.

A post-analysis isotropic reference draws512 uniform real-space spheres at the
same N. Even these have147/512,214/512,238/512 line-dominance labels
(28.7%,41.8%,46.5%). Their median L values are.311/.360/.369. Therefore a single
sparse cloud being “line-dominated” is not a significance test or a proof of a
physical filament. The dense resolved matter structure and conditional repeated
sampling are what support the interpretations above. No p-value or population
error rate is inferred from these three selected cases.

The displayed matter structure is unsmoothed beyond TSC and interpolation. At
7 Mpc/h Gaussian smoothing, mass covariance inside the same R sphere becomes
nearly isotropic (S=.925/.917/.949). This demonstrates scale dependence of this
local mass-shape statistic, not a change of tidal class. Smoothed density
covariance and the potential Hessian remain different operators. The original
native host tidal eigenvalues independently replay exactly for all three targets.

## Coordinate and provenance checks

The public CutSky X_COM/Y_COM/Z_COM origin differs by half a box from native
snapshot halo positions. Apply fixed(+1000,+1000,+1000) Mpc/h modulo2000 before
querying the matter grid. The extraction's original host_offset_mpc_h diagnostic
was the uncorrected inter-origin separation (~1732 Mpc/h), not a physical offset.
Checks on16 nearby host IDs per region reduce the median galaxy-host offsets to
.0080/.0130/.0071 Mpc/h; remaining individual offsets include HOD galaxy-host
separations. Central offsets are.0365/.0094/.0667 Mpc/h. The native host T-web
cells reproduce the stored central eigenvalues exactly, with no fitted label map.

A local affine native-to-Planck18-observer transform is estimated using positions
only and checked on234/209/312 disjoint held-out galaxies. Maximum errors are
.000723/.000746/.000778 Mpc; all below.002 Mpc. It accounts for the local observer
orientation and fiducial-distance scaling; it is not asserted as a global cutsky
transform. No truth labels enter its fit. Matter is evaluated around the galaxy
centre; the truth label remains the native host label.

Read only97^3 padded native grid cells per region, sufficient for the Gaussian
filter's4-sigma support, and interpolate to97^3 observer-frame display/quadrature
points spanning±1.5R. Large9GB raw CutSky and33GiB density sources are recorded
by path,size,mtime and build manifest; they were not fully rehashed. Small stage
manifests and saved extracted arrays are SHA256-bound. Matter below native grid
resolution is not resolved. No CFS writes, HPSS retrieval or new data transfer.

Raw CutSky redshift coordinates and exact observed coordinates agree within
9.1e-13 Mpc locally, so this paired redshift change is not a separate measurement
noise effect. Fixed observed membership replays previous N and L/P/S to1e-10.
Independent weighted-SVD validation reproduces the covariance eigensolver results,
every stored galaxy thinning draw, mass/excess shapes, probability-free counts,
and sum-to-one constraints: **4,180 checks pass**. Uniform-grid isotropy tolerance
is1e-8 (a first1e-12 tolerance failed at8e-10 due to spherical-grid discretization).
All final PNG/PDF figures visually reviewed, with matching axes and no clipped
samples. Bands are central95% thinning ranges, not confidence intervals on cosmic
population parameters.

## Execution and artifacts

CPU allocation59905073, nid200071, desi/interactive64logicalCPUs, approved30-minute
limit and Scratch. Released after12min. An initial syntax failure was corrected
before data processing. Extraction52s, analysis9s; v2 repeat7s corrected a summary
field name from dominant_fraction to dominant_counts without changing numbers.
First validator failed the over-tight isotropy tolerance; final validation/figures
passed. All successful science steps0:0; unrelated jobs untouched.

- Scripts: `workflows/p12a_vac/filament_geometry/realspace_extract.py`,
  `realspace_compare.py`, `plot_validate_realspace.py`, `realspace_sampling_null.py`.
- Evidence: `docs/evidence/filament_realspace_20261017/` (including extracted arrays).
- Figures: `docs/figures/filament_realspace_20261017/*_realspace.{png,pdf}`.
- Scratch final: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_realspace_20261017_v2/`.
- Earlier v1 is retained as a draft/provenance record; v2 is authoritative.

No new P12-A inference or model changes. No sealed phase accessed. These
controlled comparisons explain these three examples, not every wall in ph006 or
DESI. They support resolved substructure within a tidal-wall environment while
showing that sparse geometry and redshift-aperture selection can alter its label.
