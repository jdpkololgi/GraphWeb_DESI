# What the upstream BGS generator changes — 2026-09-28

## Decision boundary

Canonical v1 has enough intrinsic high-redshift galaxies to remove the old
parent-support deficit. It is NOT yet demonstrated to match observed Loa after
fibre assignment and redshift-quality selection. Do not call these scripts a
validated VAC-training replacement on the basis of raw counts.

The existing full ph000 common-sky census and the same Loa histogram give:

| z | v0.1 raw/Loa SGC | v0.1 raw/Loa NGC | v1 raw/Loa SGC | v1 raw/Loa NGC |
|---|---:|---:|---:|---:|
| .15–.25 | 1.478 | 1.361 | 1.626 | 1.466 |
| .25–.35 | 1.355 | 1.265 | 1.585 | 1.427 |
| .35–.45 | 1.135 | 1.105 | 1.612 | 1.498 |
| .45–.55 | .707 | .731 | 1.928 | 1.938 |

These are intrinsic parents vs observed data, not a matched observation model.
Both parents use numerical12<=r<19.5, no IN_Y flag, identical NSIDE256 support;
this is not a declaration that photometric selection is physically equivalent.
Sources: alignment_execution_20260926/parent/{v0.1,v1}.json and
p12a_joint_flags_20260926/HISTOGRAMS.npz. The required high-shell net retention
would be roughly52%, but that is a diagnostic ratio, NOT a thinning prescription.

## Source differences

Alex/Jade abacus branches share14222dcf794b3bff68987bf8e271203db20f8730.
The old shared_code generator and new hodpy/abacus are different prescriptions:

1. New BGS HOD fit tables replace the older population model. The default is
   `hod_fits_c000_ph000_wsys.txt`; nowsys also exists. A matching photometric
   fingerprint does not identify which HOD generated canonical v1.
2. The new LF table is predicted from the fitted HOD, with a source comment
   explicitly preserving clustering at z=.2. Its HOD defaults to
   `redshift_evolution=False`. When the E-corrected BGS target LF is used for
   evolution, its P,Q are zero. This is not equivalent to arbitrarily setting
   the old GAMA/SDSS P=1.8,Q=.7 to zero: the magnitude convention and LF differ.
3. Observed magnitudes use DESI N/S colour-dependent K corrections and E Q=.67;
   the old delivered mapping has the measured -.8(z-.1) term and GAMA tables.
   Canonical v1's DESI photometric recipe was independently replayed to near
   floating-point colour precision in the preceding report.
4. Colours have DESI-specific distributions, and satellite radial extent is
   1.5R200 in the new configuration. These can change selected clustering and
   galaxy-matter relationships; they must be tested before environmental inference.

Thus the improvement is not established as an extra IN_Y flag or an n(z)
weight. It involves both different population modelling and observation mapping.

The bundled wsys/nowsys cumulative-LF ratios at M=-20,-21,-22,-23 are
1.034,1.043,.991,1.140. That table comparison alone does not select a fit or
explain the old-to-new factor2.7 high-tail increase.

## Fixed-row diagnostic

`alignment_photometry_swap.py` samples every100th raw ph000 row in each version,
uses the same common sky and bins, and switches only the additive K/E mapping
at fixed delivered redshift/rest colour. New mapping uses the recovered
DEC32.375 split. It reports both directions. Counts are sample counts; multiply
by100 only for an approximate full-catalogue comparison. This deterministic
sample is exploratory, not an independent simulation or iid error estimate.

This does not regenerate galaxies, swap colour distributions, establish matching
absolute-magnitude conventions, or recover galaxies absent from a delivered
flux-limited parent. Consequently it measures sensitivity, not an exact causal
partition into photometry and HOD contributions. No modified catalogue is saved.

## Executed results

CPU59008849: common-sky samples166,708old /194,814v1. In .45<=z<.55:

| Fixed parent | Mapping | SGC sample count | NGC sample count |
|---|---|---:|---:|
| old v0.1 | stored old | 264 | 613 |
| old v0.1 | new K/E mapping | 750 | 1544 |
| new v1 | stored new | 822 | 1734 |
| new v1 | old K/E mapping | 314 | 704 |

At fixed old rows, changing only the additive magnitude mapping increases the
high-tail selection by2.84/2.52. Applying the old mapping to v1 removes62%/59%
of its selected high-tail rows. These changes are of the same order as the
full-census old-to-new2.73/2.65 gain. On old sampled common-sky rows, median
new-minus-old magnitude changes are -.141,-.133,-.162,-.234mag across four
shells. On v1 rows they are -.146,-.141,-.154,-.167mag (different colour mix).
A modest brightening selects many more objects in a steep bright tail.

This is direct evidence that the observation/photometry mapping can explain a
large count change without adding galaxies or changing halo occupation. It
strengthens that lead over an unsupported assertion that HOD alone is the cause.
It is not an exact decomposition: rest-colour populations, stored-M convention,
parent truncation and deterministic sample variation remain. The sparse high-z
SGC baseline differs~9% from full-census/100, so do not quote these ratios as
precision estimates or infer a fraction of the discrepancy explained.

## Requested quality/magnitude/systematics test

`alignment_v1_observation_screen.py` uses the measured old ph006 raw-to-Loa
ZWARN0 mock retention in each cap and dz=.01, including the delivered assignment
and full-HPmap selection, and transfers it to the full v1 ph000 raw histogram.
Compared against actual Loa selected counts, high-shell ratios are:

| Loa selection | SGC transferred v1/Loa | NGC transferred v1/Loa |
|---|---:|---:|
| ZWARN0, DCHI>=25, GALAXY, all targeted | 1.428 | 1.443 |
| same, 12<=r<19.5 | 1.428 | 1.524 |
| ZWARN0, DCHI>40, all targeted | 1.420 | 1.435 |
| same, 12<=r<19.5 | 1.420 | 1.515 |

The DCHI>40 alternative does NOT include the GALAXY constraint and is not
strictly a nested tightening of the baseline. All-targeted retains actual
historical targeting, including the North magnitude limit; the uniform19.5
row is only a controlled diagnostic. Historical target-bit replay already
passed for all5,436,413 selected Loa galaxies (selection-closure report).

**This transferred response fails to establish agreement.** It predicts a
43–52% high-tail excess for baseline variants. It is not an actual v1 assignment
run: denser/differently clustered targets can have different fibre losses.
Mock ZWARN0 is not a measured v1 DELTACHI2/SPECTYPE success model. Therefore
these tests cannot certify that all realistic losses resolve the discrepancy.

Explicit untested requirements for a full forward-model claim:
- v1 historical BRIGHT target replay including fibre flux, Gaia/SGA, imaging
  masks and morphology (raw v1 lacks the required observables);
- actual Loa tile/assignment and bad-fibre/exposure/HPmap mask replay on v1;
- population-dependent redshift success and quality observables, absent from
  raw v1; it would be false to claim these missing cuts passed;
- clustering/galaxy-matter validation and additional independent mock evidence.

No scalar retention is adopted as a correction. The new recipe is a supported
candidate with adequate parent support, NOT a confirmed n(z) solution. We can
continue independently by constructing a pinned observational forward model;
a producer reply is not a prerequisite, but missing observables cannot be
replaced by an undocumented pass flag. No fitting, regeneration or VAC change.
