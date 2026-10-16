# Is tracer density the main cause? Abacus follow-up

Date: 2026-10-15 (client date). Follow-up to
`GraphWeb_DESI/docs/legacy_gat_sequential_tests_20261014.md`.

**Conclusion:** sparse tracers are a demonstrated information limitation in
both simulations. They have not been established as the main cause of the
old/current observed property-contrast difference. Simulation mass resolution,
galaxy population, observational selection and target smoothing must be varied
separately. The earlier suggestion that Abacus is simply less dense because of
its lower mass resolution was not supported by the density comparison.

## Fresh Abacus measurement

Reused the exposed historical ph000 S2 observed/harmonized wedge, not reserved
phases or current DESI predictions. Truth is the sum of saved R7 eigenvalues.
Counts are self-excluded, in observer-coordinate apertures; the query galaxies
stay fixed under thinning. A 14 Mpc/h buffer excludes the rectangular angular
and radial boundaries. This does not certify internal survey-mask completeness.
The first three shells use 10,000 deterministic queries each; the last uses
all 2,424 eligible queries. Counts use the final selected shell population,
not the denser unharmonized buffer catalogue.

| Redshift | rho(N2, true trace R7) | rho(N7, true trace R7) | rho(selection-adjusted N7, trace R7) | Median N7 | No N7 neighbours |
|---|---:|---:|---:|---:|---:|
| .15–.25 | .341 | .760 | .791 | 22 | .34% |
| .25–.35 | .247 | .708 | .754 | 8 | 3.36% |
| .35–.45 | .169 | .620 | .677 | 2 | 17.40% |
| .45–.55 | .077 | .414 | .496 | 0 | 61.14% |

rho denotes Spearman rank correlation. Selection adjustment is
log[(N+.001)/(ntilde(z)*V+.001)], using the frozen historical spline. Pooled
shell correlations can mix radial trends with spatial information; the receipt
also contains .02-wide redshift slices. For example, raw N7 correlations in
the first shell's narrow slices range .758–.831. The all-row historical N7
correlations are .732/.694/.594/.397; differences from the primary table reflect
interior/sample selection, not a changed field.

**These numbers are not a matched comparison to TNG's .891 or .723–.733:** the
TNG experiment used periodic real-space tracers and native R~2 truth; this
Abacus measurement uses an observed wedge and R7 truth. In particular, N2 in
the table is a 2 Mpc/h count aperture compared with **R7** truth. It is not a
new test of R2 truth. Also, trace is the density component of the tensor;
high count/trace correlation does not establish accurate tidal shear or full
web classification. The corresponding all-row N7/lambda1 correlations are
.550/.533/.492/.347, substantially lower than those for trace.

Randomly thinning the same Abacus tracer population gives:

| Redshift | Full density: rho(N7, trace R7) | 30% retained | 10% retained |
|---|---:|---:|---:|
| .15–.25 | .760 | .688–.703 | .566–.579 |
| .25–.35 | .708 | .585–.601 | .423–.434 |
| .35–.45 | .620 | .433–.458 | .289–.302 |
| .45–.55 | .414 | .231–.256 | .120–.137 |

Ranges are extrema over three random seeds, not confidence intervals. These
are within-shell interventions at fixed truth and queries: they identify an
effect of random tracer removal on this statistic. They do not quantify the
effect of luminosity selection, different HOD bias, changed matter resolution,
or the observed old/current property gap.

Validation: fresh N7 counts exactly reproduce cached integer counts on all
32,424 queries; 48 independent brute-force aperture checks pass. The first
source inspection found 50 negative BOX_INDEX rows in the lowest nominally
valid shell, so their truth labels are excluded explicitly. They remain
observed tracers, since counts require positions rather than their labels.
The other shells have no such exclusions. All input files and the script are
hashed in RESULTS.json. No training, catalogue change, graph features or
reserved-phase access occurred. Slurm inspection remains unavailable; this
short one-thread diagnostic used the previously authorized login fallback.

Reproducer: `workflows/abacus_tweb/count_trace_density_probe_20261015.py`.
Results: `/global/u2/d/dkololgi/fmsu_tokyo_20261010/abacus_density_followup_20261015/RESULTS.json`.

## Existing stronger intervention: adding real Faint neighbours

The saved MT4 closeout, freshly inspected in
`docs/evidence/p8/multitracer_mt4_decision.json`, is more informative than a
comparison of two unrelated catalogues. It keeps Bright target galaxies and
truth, adds Faint context, and compares it with angularly scrambled Faint
context retaining the nuisance sampling/selection structure.

| Model/context | Macro R2(lambda1), rotation 0 | Rotation 2 |
|---|---:|---:|
| Bright only | .5070 | .5197 |
| Bright + actual Faint positions | .6074 | .5394 |
| Bright + scrambled Faint positions | .4753 | .5187 |

Actual minus scrambled Faint improves every shell in both rotations. Actual
minus Bright improves the sparse shell by .1392/.1091. This establishes
additional spatial information within ph000, beyond merely increasing input
occupancy. It is not proof that generic random density increases are equivalent
to adding Faint galaxies, which are a different tracer population. The model
failed the rotation-2 promotion safeguards, and there is no new independent-
phase or old-versus-current property-gap conclusion.

## Why R7 was chosen

The user clarifies that R2 was tried and its noisy/weakly recoverable field
motivated R7 at lower tracer density. Treat R7 as a deliberate information–scale
tradeoff, not an accidental mismatch or evidence of poor model design. The
June 26 science record independently supports this framing: scale-matched
apertures improve recovery, while larger smoothing erases compact clusters.
The precise original R2 run receipt was not recovered in this follow-up.

Two forms of noise must remain distinct. A sparse **galaxy-derived estimator**
can be noisy while the true matter field is numerically converged. If the
**matter-derived truth itself** is noisy at R2, test particle subsampling,
mesh/deposition and interpolation convergence separately. Thinning galaxies
does not change the fixed true matter field. The reason to include R2/R7 in
the next experiment is to separate a legitimate adaptation from its scientific
consequence: different scales describe different environments and can have
different true property associations even with perfect inference.

## Mass resolution, HOD and the observed tracer population

Mass resolution sets which host haloes and their structure are explicitly
resolved; stellar mass is a different quantity. HOD sampling then assigns
central/satellite galaxies conditional on hosts. The observational selection
determines which of those galaxies are retained as usable tracers.

The published Abacus BGS construction uses a resolved-host cut of
10^11 Msun/h and supplements lower-mass hosts with field-particle proxies.
Its cubic Mr<-18 population has density .037 h^3 Mpc^-3, already comparable
to this archived TNG graph's .03296. Those are published-parent versus local-
TNG values, not a fresh census of our local Abacus parent.
[Smith et al., sections 3.3 and 5](https://arxiv.org/html/2312.08792v2#S3.SS3).

Thus HOD modelling can supply faint galaxies statistically, including satellites
without requiring one resolved subhalo per satellite. It cannot make galaxies
below the survey flux limit observable, nor guarantee that unresolved-host
placements reproduce every galaxy–matter or environmental statistic. More
particles can improve that placement model without increasing the retained
BGS density when the luminosity function and survey selection are held fixed.

The important diagnostic is the full chain: resolved matter/hosts → populated
galaxies → magnitude-selected targets → observed tracers. Compare densities
and host/proxy fractions at each stage. The previous 7.6 ratio compared a TNG
mass-selected box with an observed-survey selection-density estimate; it did
not isolate any one stage of that chain.

## How to establish whether density is the main cause

Begin with the current artifacts, then use controlled experiments rather than
another architecture ranking. Do not use the disjoint observed samples to
estimate a causal percentage.

1. **Fix the physical target.** In one high-resolution box at one epoch, produce
   independently converged R2 and R7 truths from the same matter field. Keep
   query galaxies, galaxy properties and spatial validation regions fixed.
   Measure the true property contrast at each scale before learning anything.
   This distinguishes genuine smoothing of the property relation from inference
   attenuation.
2. **Intervene only on density.** For each truth scale, use dense, intermediate
   and survey-density random tracer subsets, rebuilding inputs each time. Train
   the same model with matched budgets/seeds, and also evaluate a fixed model
   under thinning. The first measures achievable recovery; the second measures
   sensitivity to distribution shift. Score eigenvalues, trace, shear, classes,
   calibration and predicted versus true property contrast on identical queries.
3. **Add observational effects one at a time at fixed density.** Compare random
   thinning with luminosity selection, then RSD, mask and assignment effects.
   Equal number density does not imply equal galaxy bias or spatial information.
   This tests whether selection, rather than sparsity alone, explains a residual.
4. **Change simulation only after matching the experiment.** Repeat at common
   epoch, target, density and selection in verified Abacus products. For mass
   resolution specifically, use the matched high/highbase phase pair with
   harmonized HOD/selection and number density; do not compare unrelated TNG
   and Abacus volumes and call the difference resolution.

A density-dominated explanation would require thinning alone to reproduce
most of the lost inference skill and attenuation of predicted property
contrasts within this controlled design, with small residual changes from
selection and simulation choice. Quantify the fraction for a predeclared
error/attenuation measure with paired spatial uncertainty; do not subtract
Spearman coefficients from different targets or divide the observed 21.87/8.97
pp contrast gap into a causal budget. If true R7 property contrasts are already
weak before inference, smoothing/estimand choice contributes independently.
If equally sparse samples differ strongly, galaxy bias, HOD or observation
modelling remains material. These neural/paired-field tests are not yet run.
