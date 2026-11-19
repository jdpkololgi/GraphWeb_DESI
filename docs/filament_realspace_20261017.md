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

## Launch contract

Run these scripts with the working directory
`workflows/p12a_vac/filament_geometry`. `realspace_extract.py` imports
`prepare_parent` and `mock_truth` by bare filename. The other three also
import `compare_offslab` that way. There is no package.
`mock_truth.ROOT` is `/pscratch/sd/d/dkololgi/abacus/p10_multiphase` and
`mock_truth.R` is `7/0.6766` Mpc. Observer coordinates follow
`prepare_parent.coordinates`: Planck18 comoving distance, with the RA origin
at 180°. This is the exposed ph006 development mock. The completed v2 tree is
the scratch path below; a new run needs a new directory.

```bash
cd workflows/p12a_vac/filament_geometry
python realspace_extract.py --out /pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_realspace_NEW
python realspace_compare.py /pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_realspace_NEW
python realspace_sampling_null.py /pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_realspace_NEW
python plot_validate_realspace.py /pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_realspace_NEW
```

`realspace_extract.py` takes `--out` and raises `FileExistsError` if that
path exists. It checks `SLURM_JOB_ID` only; it does not require an `nid`
hostname. Targets are the `targetid` list in
`docs/evidence/filament_mock_examples_20261017/EXAMPLES.json`, in that order.
Each centre is FITS row `TARGETID−1` of the bright parent
`ph006_bgs_bright_parent_with_tweb_eigs_rs7_ngrid2048_thr0p2_15d.fits`, and
the read is asserted to return that `TARGETID`. The extraction keeps raw
CutSky and parent rows within 45 Mpc in real space or in redshift space.
Analysis still uses `R`. The public raw file is
`cutsky_BGS_z0.200_AbacusSummit_base_c000_ph006.fits`. Large inputs are
recorded by path, size, and mtime with `full_hash_checked` false. The
per-region npz and `EXTRACTION.json` land in `--out`.
`host_offset_mpc_h` in that JSON is the uncorrected inter-origin separation.
The `+1000` inside that expression is not the snapshot-origin shift applied
later in `realspace_compare.py`. A pure L/2 offset still has length about
1732 Mpc/h there. The local affine map is fit on every other galaxy that is
inside 30 Mpc in both frames and more than 1 Mpc from the centre in the
observer frame. The held-out half is the check. The central host read needs `abacusnbody` and that
galaxy's `halo_info_*.asdf`.

`realspace_compare.py` takes the extraction directory as a positional
argument. It requires `SLURM_JOB_ID`, reads `EXTRACTION.json`, and raises
`FileExistsError` if `COMPARISON.json` is already there, before the density
cube is opened. Region npz files are replaced on a rerun that has not yet
written `COMPARISON.json`. Before importing
`workflows.abacus_tweb.p12a_coordinate_sample_audit` it inserts the hardcoded
path `/global/u2/d/dkololgi/TNG/Illustris` at the front of `sys.path`. It
does not read `ILLUSTRIS_ROOT`. The density cube comes from the first
`ph006/targets/density/*.manifest.json`; the cube must be `2048³` with
`processed_file_count` 136. One `default_rng(171019)` is created before the
region loop, so galaxy thinning and the later matter draws share one stream
in `EXTRACTION.json` order. Thinning is 512 draws without replacement and is
marked unavailable when the real-space pool is smaller than the observed
count. Matter draws are 512 samples with replacement on the display grid
inside `R`, weighted by the interpolated TSC counts. The smoothed mass field uses
`gaussian_filter` with σ = 7/Δx cells, `truncate=4`, and `mode='nearest'`
on the 97-cell patch. The script asserts a held-out transform residual below
0.002 Mpc, a median origin-corrected host offset below 1 Mpc/h, an exact
replay of the stored central eigenvalues, and the previous redshift-space
shape to `1e-10`. Exact observed positions are the canonical `points.npy`
rows reordered as `(-y, -x, z)`.

`realspace_sampling_null.py` has no `__main__` guard. Importing it, or
running it, reads `sys.argv[1]`, requires `SLURM_JOB_ID`, and overwrites
`ISOTROPIC_NULL.json`. Seed `171020` is one stream shared across regions.
Each draw is a uniform ball of the observed count. The shape components are
asserted to sum to 1 at `1e-12`.

`plot_validate_realspace.py` does not check `SLURM_JOB_ID`. It overwrites
`VALIDATION.json` and `{key}_realspace.{png,pdf}` in the results directory.
The copies under `docs/figures/filament_realspace_20261017/` are the archived
figures from the completed run. The cross-check uses an independent weighted
SVD: galaxy shapes at `1e-12`, mass shapes at `2e-7`, and the empty-sphere
isotropy residual at `1e-8`. The first validator failed a `1e-12` isotropy
tolerance at about `8e-10` because of the spherical-grid discretization. The
script updates the in-memory report with the raw-versus-observed redshift
residual and does not rewrite `COMPARISON.json`.

## Execution and artifacts

CPU allocation 59905073, nid200071, desi interactive, 64 logical CPUs, approved
30-minute limit, and Scratch. Released after 12 min. An initial syntax failure
was corrected before data processing. Extraction 52 s, analysis 9 s; the v2
repeat (7 s) corrected a summary field name from `dominant_fraction` to
`dominant_counts` without changing numbers. The first validator failed the
over-tight isotropy tolerance; final validation and figures passed. All
successful science steps exited 0. Unrelated jobs were left untouched.

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
