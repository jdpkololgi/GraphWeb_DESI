# Environment-trend investigation: tests, simulation choices and transfer learning

Date: 2026-10-14. This is the consolidated report requested for the old
IllustrisTNG GAT versus current Abacus-trained P12-A investigation. “Completed”
means the stated diagnostic was executed; it does not imply its hypothesis
was established. “Open” identifies work for which there is no result yet.

The completed checks isolate a substantial dependence on stellar-mass mixture
and probability weighting. Switching CG15 to CG5 preserves the P12-A trends.
The tested mass-completeness cut barely changes the P12-A contrast after
holding the control mixture fixed. The old label field is consistent with
about 2 Mpc/h smoothing, while P12-A targets 7 Mpc/h. A remaining difference
cannot yet be assigned to model quality: populations, physical targets and
training domains differ. A new controlled TNG thinning experiment supports
tracer density as an information limitation, but does not measure its causal
contribution to the observed property-contrast gap.

The evidence favours retaining Abacus for independent-phase survey training
while adding a high-resolution comparison/pretraining branch. No researched
alternative is currently verified to provide all three of: many independent
high-resolution phases, suitable BGS survey mocks, and a complete galaxy-to-
matter truth mapping for our full redshift range. The strongest immediate
alternative comparison is Uchuu; Rosella is especially relevant to BGS
galaxy–environment modelling. Neither is a demonstrated drop-in replacement.

## Test and findings ledger

Unless specified otherwise, contrasts are cluster minus void, “soft” means
probability weighting, and quenching means log10(sSFR/yr^-1)<-11. Percentage
points are abbreviated pp. Different rows can have different common support;
their effects must not be added as a causal budget.

| Test | Status and sample | Finding | Interpretation / remaining limit |
|---|---|---|---|
| Saved old/current catalogue overlap | Completed; all IDs | Seven shared IDs, zero consistent redshifts; old z=.01–.06, current .15–.55 | Existing VACs cannot support a same-galaxy model comparison. |
| Old duplicate-ID audit and CG15 join | Completed; 199,314 old rows | 487 repeated IDs excluded from controlled comparisons; 193,706 unique, consistent CG15 matches | Avoid treating conflicting duplicate observations as interchangeable. |
| Old north/south ordering | Completed; all 199,314 rows | IDs, positions, redshifts and masses exactly reproduce corrected north-then-south order; historical fix 5296648 | No recurrence of that bug in this saved VAC. |
| Current inference-to-VAC ordering | Completed; 5,436,413-ID census, 1,009 rows in eight cores, 256 property-cache rows | Exact identities, coordinates, saved-draw summaries and merge parity; 255 fresh CG15 source matches | Strong checks against row reassignment; not a full neural-network replay. |
| Common CG15 estimator and mass/z control | Completed; 187,358 old / 75,321 current controlled rows | Soft quenched contrasts 21.87 / 8.97 pp | The old association remains stronger on each sample's own controlled support. |
| Shared between-sample mass reference | Completed; 63,373 / 74,816 rows | Soft 14.37 / 8.81 pp; hard 14.89 / 11.85 pp | Much of the original apparent gap depends on mass mixture and probability weighting; redshift remains unmatched. |
| Original completeness rationale | Completed reconstruction; parent before mass cut | Faintest-20% / 90th-percentile limit reaches logM=9.000 at .055–.06 | Supports the remembered reason for the low-z ceiling; exact historical decision receipt missing. Raw FSF>=9 corresponds to physical logM>=9.305 with h=.704. |
| Main-BRIGHT restriction | Completed | Old soft contrast 21.872→21.925 pp | Main versus ancillary target mix is small in this test. |
| Common angular wedge | Attempted, inadequate controlled support | 7,305 eligible old galaxies, only 172 survive strict shared class/mass/z cells | No robust sky-only attribution. |
| Current 90% limiting-mass cut | Completed; 52,243 parent / 49,403 selected on common cells | Soft 9.459→9.425 pp; hard 11.852→12.041 pp | Does not recover a large missing contrast at fixed control mixture. |
| Stronger completeness limit +.2 dex | Completed; 24,849 / 22,951 | Soft 8.940→8.562 pp; hard 10.251→10.146 pp | Same conclusion for this sensitivity cut; other selection effects remain. |
| Probability sharpening on identical rows | Completed | Old p→hard 21.872→22.667 pp; current 8.970→11.764 pp | Pure estimator sensitivity, not a calibration repair. |
| Confidence selection | Completed at maxP>=.7; higher current thresholds fail support | Current paired soft 9.735→11.119 pp; old 22.014→23.455 pp | Stronger selected associations; maxP>=.85/.95 cannot be compared reliably for current four classes. |
| CG5 versus CG15 | Completed; 75,321 identical controlled rows | All seven trend directions persist; soft quenching 8.970→9.573 pp, log-sSFR -.2274→-.2177 dex | SED model choice is not the main explanation here. |
| Identify the old checkpoint | Completed; eight queries, eight candidate files | January 15, 2026 checkpoint reproduces probabilities to 1.19e-7; four/five-hop parity | Identifies this saved VAC, not necessarily the September 2025 slide. |
| Recover old target definition | Completed field transfer and label probes | Effective Gaussian R=2.02662 Mpc/h; historical R=2; threshold .2; 8,192/8,192 sampled labels match | Old and P12-A R=7 labels are different physical targets. |
| Training-cache versus field labels | Completed; both 283,969-node caches | 283,967 labels agree; two discrepancies remain | Very strong scale/provenance evidence, not perfect voxel-level closure. |
| Independent density versus current classes | Completed earlier; 91,066 common rows | Soft 7.42 pp; N4 7.82; N7 7.68; projected fifth neighbour 9.52 | Current trends agree in sign and scale with non-neural density estimates. |
| Additional N7 control | Completed earlier; 4,539 common rows | Soft contrast 4.56→2.17 pp | Density accounts for substantial association; not a causal mediation result. |
| Independent density in old wedge | Executed; 2,941 buffered queries | Counts verified, but shared class/mass/z support fails | No controlled old/current density-amplitude conclusion. |
| TNG random thinning at fixed truth | Completed; 10,000 fixed queries, three seeds | At 10% density, R2 count/trace Spearman falls .891→.723–.733; zero-neighbour fraction 7.6%→45.8–46.3% | Direct evidence for loss of environmental information with sparsity; no retraining or property outcomes in this test. |
| Low-z product and epoch inventory | Completed bounded inventory and 512,000-row probe | Low-z BOX_INDEX=-1 occurs in resolved and unresolved rows; inspected cutsky products all labelled snapshot .2; base archive ends at .1 | BOX_INDEX=-1 does not identify a different epoch or only unresolved hosts. |
| Low-z raw halo-index mass match | Completed; one slab, 12,000 low-z rows | Snapshot .2 matches 0.37%; .1 matches 0.32%; high-z .2 control matches 92.7% | Changing snapshot alone fails. |
| Reconstruct filtered low-z halo index | Executed, failed; same slab | Published mass/radius cuts match only .53% raw / .55% cleaned masses | Public template explains different index semantics but is not an exact production crosswalk. |
| Same-row, same-truth GAT/P12-A comparison | Open | Checkpoint recovered, but compatible graph/normalization and paired R2/R7 truth not assembled | No architecture winner or model-causal explanation claimed. |
| Alternative simulation audit | Completed literature plus existing local product records | See candidate table below | Phase count, matter access and BGS mock maturity assessed separately. |
| High-resolution pretraining then Abacus adaptation | Designed, not trained | Factorial comparison and acceptance criteria below | Promising hypothesis, no demonstrated gain yet. |

Results live in
`/global/homes/d/dkololgi/fmsu_tokyo_20261010/legacy_selection_tests_20261014/`.
CG5 results and figure are in the neighbouring `p12a_cg5_check_20261014/`.
Initial diagnostics used the authorized 15-minute, single-thread login batch.
The subsequent user instruction authorized additional compute and a login
fallback. An interactive CPU allocation was attempted but failed before job
creation because this tool environment could not establish Slurm configuration.
Subsequent diagnostics used one CPU thread with bounded memory. No training,
full mock inference or production-catalogue modification was performed.

## 1. First controls, completed earlier

The same CG15 source and plot recipe give old/current probability-weighted
cluster-minus-void quenched-fraction contrasts of 21.87/8.97 percentage points.
Matching the between-sample stellar-mass distribution gives 14.37/8.81 pp;
hard labels on that common-mass comparison give 14.89/11.85 pp. The apparent
gap shrinks from 12.90 to 3.04 pp across these choices. These stages change
support and weights; the reduction is not an additive causal decomposition.
The samples have disjoint redshifts, so this is not a same-galaxy comparison.

## 2. Selection and completeness

Reconstructed the actual main BGS_BRIGHT targeting limits: extinction-corrected
r<19.5 South and r<19.54 North for current Loa. All 103,788 targeting-cache
IDs, sky positions and redshifts match the current property cache exactly.
Old comparisons use main BRIGHT objects in the same angular wedge.

In each redshift slice, take the faintest 20% relative to the targeting limit,
then the 90th percentile of logM+0.4(r-r_limit). The old curve uses the parent
catalogue before its mass cut, FastSpecFit masses converted with h=.704, and
dz=.005. It reaches logM=9.000 in .055<z<.06. All 7,305 eligible old main-BRIGHT
wedge galaxies pass this curve and a further +.2-dex threshold. This supports
the remembered low-z completeness rationale; it does not recover the exact
original numerical selection receipt. The old loader's raw LOGMSTAR>=9 is
log physical mass>=9.305 under that conversion, rather than exactly 1e9 Msun.

For current CG15, the limit rises from logM=10.496 at .20<z<.21 to 10.884 at
.29<z<.30. The cut keeps 59,284 of 102,622 eligible galaxies; adding .2 dex
keeps 33,174. Outcome/control masses are CG15 in both samples. Within each
parent-versus-cut test, the mass/z cells and reference mixture are identical,
with >=20 effective galaxies in every hard and probability-weighted class.

| P12-A paired test | Parent soft contrast | Selected soft | Parent hard | Selected hard |
|---|---:|---:|---:|---:|
| 90% limiting-mass curve | 9.459 pp | 9.425 pp | 11.852 pp | 12.041 pp |
| Curve plus .2 dex | 8.940 pp | 8.562 pp | 10.251 pp | 10.146 pp |

The first test retains 52,243 parent / 49,403 selected galaxies on common
support; the stronger cut retains 24,849 / 22,951. These cuts do not restore
the much larger old contrast. They diagnose flux-limit mass selection, not
surface-brightness selection, fibre assignment, redshift failures or SED bias.
Predictions and the tracer population used to infer them remain fixed.

Old main-BRIGHT versus all-old contrasts barely change (21.872 to21.925 pp).
Matching its angular wedge leaves only172 selected galaxies on the strict
common cells: its apparent contrast change is not a robust sky-effect estimate.
Do not interpret that small-support comparison as successful sky matching.

## 3. Probability weights and confidence selection

On exactly the same rows/control mixture, sharpening probabilities changes
the P12-A contrast from8.970 pp to10.695 (normalized p-squared),11.468
(normalized p-fourth), and11.764 (hard labels). Old GAT changes from21.872 to
22.511,22.672,22.667 pp. Thus membership weighting alone changes P12-A's
reported association by2.79 pp and old GAT's by0.80 pp. Sharpening is an
arithmetic sensitivity check, not a calibration correction.

For maxP>=.7, P12-A's paired parent/subset soft contrasts are9.735/11.119 pp
(hard12.243/14.392), retaining43,189 parent and27,561 selected rows. Old GAT
gives22.014/23.455 pp (hard22.815/24.241). P12-A maxP>=.85 and.95 fail common
four-class mass/z support. Old high-confidence subsets have stronger
contrasts, but these are selected populations, not evidence of greater accuracy.

CG5 independently preserves all seven P12-A property-trend directions on
75,321 identical controlled galaxies. Soft quenched contrast is8.970 pp with
CG15 versus9.573 pp withCG5; log-sSFR contrasts are-.2274/-.2177 dex. The
paired CG5-minus-CG15 contrast change is+.603 pp (16–84% sky interval
.437 to.793 pp). Alternative CG5-mass controls also preserve the result.

## 4. Physical target scale, checkpoint and independent density

Recovered eight archived GAT checkpoint candidates. The January15,2026 file
`trained_gat_model_ddp_2026-01-15.pth` reproduces all eight sampled old-VAC
probability vectors to max error1.19e-7. Exact four-hop inference used4,325
nodes; adding a fifth hop changes results by at most1.19e-7. Catalogue IDs
and redshifts agree exactly. This identifies the checkpoint behind the saved
VAC, not necessarily the original September2025 slide.

The saved TNG density-to-tidal-trace Fourier transfer has effective Gaussian
radius2.02662 in the 205-Mpc/h box, consistently on all three axes, with
maximum relative fit error2e-5. Historical commit5fd247b specifies Rsmooth=2
for those saved field paths; the slight effective-radius offset is unresolved
numerical smoothing, not evidence of R=7. Threshold0.2 reproduces all8,192
sampled field classifications. Both archived training graph caches reproduce
283,967 of283,969 field labels when their positions are mapped back to this
grid; two voxel-assignment discrepancies remain. Cache target arrays agree
with graph labels exactly. This is strong numerical and code evidence for
the old approximately2-Mpc/h label definition, with the two-row caveat retained.

The ten old model inputs are degree, three edge-length summaries, clustering,
density, neighbour density, and three inertia eigenvalues. There is no direct
stellar-mass or colour input. Cached transformed features were used for replay;
reproducing that normalization on a new population remains a transfer issue.
The training code splits nodes in one graph; this does not establish performance
on a held-out survey volume.

The existing current-P12-A density experiment already measures quenched
contrasts7.42 pp for soft classes,7.82/7.68 pp for4/7-Mpc/h counts, and9.52 pp
for projected fifth-neighbour density on the same91,066 controlled galaxies.
With an additional N7 control, the soft-class contrast changes from4.56 to
2.17 pp on the same4,539 galaxies. Much of the class association follows
measured density; this is not a proof of causal mediation or zero tidal effect.

A new old-sample1/2/4/7-Mpc/h count check used198,340 unique old tracers and
2,941 property-valid queries after a7-Mpc/h angular-wedge/radial buffer.
Twelve brute-force neighbour checks pass. It fails shared four-class support
at both .1 and .2 dex mass binning, including after removing tied N1 quartiles.
Therefore no old-versus-current density amplitude is claimed from this probe.
The rectangular buffer also does not certify the old survey mask completeness.

## 5. Genuine model comparison: prerequisites recovered, execution open

A valid matched-target model comparison has NOT been completed. Existing
observational VACs have no consistent same-galaxy overlap. The recovered GAT
and P12-A also target different smoothing scales. The latest compute approval
allows further work; authorization is not the outstanding limitation. A
scientifically interpretable comparison still needs compatible model inputs,
normalization, common spatial support and both independently verified truths.

Concrete next experiment: use one exposed Abacus validation phase (ph006),
one frozen main-BRIGHT tracer selection, identical spatial core/context and
row IDs, and independent truth at BOTH R=2 and R=7 with threshold0.2. Build
the old ten-feature alpha graph with its recovered normalization and run the
identified GAT checkpoint alongside the frozen P12-A model. First report
the full2-model by2-truth-scale confusion/calibration matrix. This measures
frozen-method transfer and scale sensitivity, not architecture alone. A fair
architecture comparison then requires both methods trained for the same truth
scale with identical phase separation. Do not choose the winner by observed
property-trend strength. Colour in the Abacus mock is HOD-assigned, so its
property contrasts are not an independent validation of galaxy-formation physics.

## Added low-redshift provenance test

The local inventory contains BGS CubicBox/v0 at z=.1,.2,.4, but v0.1 and v1
only at z=.2. CutSky/BGS v0,v0.1,v1 each list z=.2 only. The public ph000 halo
archive begins at z=.1 and contains no z=0 snapshot directory. This is a
statement about these inspected products, not all Abacus data worldwide.

The July record's claim that snapshot-shell stitching was externally verified
is too strong. [Smith et al., sections5.2 and5.4](https://arxiv.org/html/2312.08792v2#S5.SS2)
describe a single z=.2 snapshot for their BGS cutsky, with evolving magnitudes
and unresolved low-mass hosts placed on field particles. This paper does not
establish the exact provenance of every row in our internal v0.1 product.

A fresh512,000-row source probe finds BOX_INDEX=-1 below.1 for BOTH resolved
and unresolved populations; it is not solely an unresolved-host flag. Around
.15 there is a transition, with some mixing when binned by cosmological
redshift. In slab FILE_NUM=1, the halo mass associated with HALO_INDEX matches
the stored mass within1% for92.7% of6,608 z=.2–.3 controls at snapshot.2,
but only0.37% of12,000 z<.15 rows. Substituting snapshot.1 matches only0.32%
of those low-z rows. This fingerprint is a diagnostic, not a complete ID or
coordinate crosswalk; the higher-z unmatched remainder also needs explanation.

The upstream source now provides a more specific explanation. The pinned
[older generator](https://github.com/amjsmith/shared_code/blob/a40a5d3be570a6ce9d5a0698268b04864d4d3af2/abacus_mocks/make_mocks_script.py)
uses snapshot .2 for both branches: an ordinary cubic population and a
separately populated, fainter low-z population. Before populating the latter,
it filters resolved haloes by mass and observer-centred radius. Its
[galaxy catalogue implementation](https://github.com/amjsmith/shared_code/blob/a40a5d3be570a6ce9d5a0698268b04864d4d3af2/abacus_mocks/hodpy/galaxy_catalogue.py#L110-L114)
then assigns halo_ind using consecutive indices into that already filtered
table. The [newer cutsky implementation](https://github.com/amjsmith/hodpy/blob/14222dcf794b3bff68987bf8e271203db20f8730/make_mocks/cut_sky_tools.py)
explicitly assigns box_ind=-1 to the low-z branch; high-z box_ind refers to a
cubic **galaxy** row. Merging uses observed redshift .15, explaining the mixed
transition when the catalogue is instead binned by cosmological redshift.

This refutes the inference “BOX_INDEX=-1 therefore a z=0 or different-snapshot
population.” It does not certify that the public template exactly generated
our local version. Reconstructing its log halo mass>=11 and zmax=.16 radius
cuts using the c000 header cosmology, observer (-1000,-1000,-1000), and
x_L2com still matches only .53% (raw) / .55% (cleaned) of 12,000 stored masses.
The source template itself uses c004, so it is a template, not a c000 run
receipt. That negative result prevents us from declaring the mapping fixed.

The next repair must recover the exact production halo/filter/replication
mapping, or independently validate the true-position-to-periodic-box transform
for the resolved and particle-based populations. Sample the appropriate .2
matter field only after that transform is established, and rerun coordinate,
density-correlation and permutation controls. Do not substitute snapshot .1,
invent z=0 products, or lower the science cut. The observed invalid labels
remain sufficient reason to retain the current exclusion.

Receipts: RESULTS.json, PREPARATION.json, checkpoint_replay.json,
label_scale.json, training_label_parity.json, old_density_scale.json,
lowz_product_probe.json, lowz_filtered_index_probe.json, tracer_thinning.json.
Public source snapshots are retained as inert text with pinned upstream
commits; they were inspected, not executed. Source hashes are collected in
MANIFEST.json.

## Tracer density: what the new experiment establishes

The archived TNG graph has 283,969 galaxies in (205 Mpc/h)^3, a number density
of .03296 h^3 Mpc^-3. The current wedge's recorded selection-density field has
10th/50th/90th percentiles .00265/.00436/.00681 h^3 Mpc^-3 after unit conversion.
Thus TNG is about 7.56 times denser than the median of that current field.
This is a comparison of a periodic mass-selected population with a survey
selection-density estimate, not a claim of equal galaxy bias or completeness.

At the same 10,000 query galaxies, keep the native R~2 tidal trace fixed and
independently thin the TNG tracers with three random seeds. Exclude each query
from its own neighbour count, whether or not it survives thinning. Periodic
boundaries are used, and brute-force checks pass. The following ranges are
the extrema across three thinning seeds, not confidence intervals.

| Retained fraction | Approximate density (h^3 Mpc^-3) | R2 count/trace Spearman | R7 count/trace Spearman | Fraction with zero R2 neighbours |
|---:|---:|---:|---:|---:|
| 1 | .03296 | .891 | .866 | .076 |
| .3 | .0099 | .824–.830 | .848–.850 | .240–.247 |
| .1 | .0033 | .723–.733 | .805–.809 | .458–.463 |
| .03 | .00099 | .588–.597 | .699–.709 | .692–.695 |
| .01 | .00033 | .431–.447 | .547–.550 | .844–.848 |

This is direct evidence that denser tracers make the same underlying field
easier to recover. Larger apertures are more resistant to thinning, while
averaging over a different spatial scale. R7 in this table is a **count
aperture**, still compared with the native R~2 truth; no R7 truth was generated.
Random thinning does not emulate luminosity-dependent bias, fibre assignment,
redshift-space distortions or the survey mask. It therefore supports the
density mechanism without assigning a numerical share of the old/current
property-contrast difference to it.

Three further reasons prevent interpreting stronger observed trends as better
field reconstruction: old versus current targets have different smoothing;
low-z galaxies and the current flux-selected population differ; and stronger
association with a density statistic is possible even for a poorly calibrated
web posterior. The old GAT uses geometry/density features, not mass, colour or
SFR directly. A model need not learn galaxy properties explicitly to select
environments in which those properties differ.

## Simulation alternatives: evidence checked on 2026-10-14

Independent phases are distinct initial-condition realizations at the same
cosmology. Rotations, observers, repeated HOD draws, resolution variants and
feedback variants are not additional independent phases. A matched
fixed-and-inverted pair is valuable but is not two independent Gaussian draws.
Mass resolution alone does not guarantee a denser or more realistic BGS sample.

| Candidate | Volume / resolution | Independent-phase situation | BGS and matter-truth readiness | Assessment for this task |
|---|---|---|---|---|
| AbacusSummit base | 2 Gpc/h; 6912^3 particles | 25 c000 base phases | BGS mocks and existing pipeline; particle subsamples and halo products; known low-z mapping and high-z population problems | Best verified ensemble foundation. Repair the mock-to-field interface and selection model rather than discard established phase separation. [Specifications](https://abacussummit.readthedocs.io/en/latest/simulations.html), [products](https://abacussummit.readthedocs.io/en/latest/data-products.html). |
| Abacus high / highbase | 1 Gpc/h; high has 6300^3 particles and 6x finer mass resolution than base | One matched phase ph100 at high and base resolution; not an ensemble | Full snapshot outputs, no simulation lightcone; ready BGS product for this pair not verified | Particularly clean resolution-control or pretraining candidate within the same simulation family. Requires BGS population construction/verification. [Specifications](https://abacussummit.readthedocs.io/en/latest/simulations.html). |
| Abacus small | 500 Mpc/h at base mass resolution | 1,643 realizations reach z=.2 | Particle subsamples available; not a verified ready BGS training set here | Useful if phase diversity, rather than resolution, is limiting. Does not solve faint-halo resolution. [Specifications](https://abacussummit.readthedocs.io/en/latest/simulations.html). |
| Uchuu / Uchuu-DESI | 2 Gpc/h; 12800^3; particle mass 3.27e8 Msun/h | One flagship volume. The 102 public EDR lightcones are views of that volume, not 102 independent gravity simulations | Published SHAM BGS mocks, evolving snapshots and halo/subhalo structure. Galaxy-to-field crosswalk for our local products not established | Strongest alternative BGS/high-resolution comparison. Cannot replace independent-phase validation alone. Local Y3-v2.0 products previously tested have no .50–.55 support. [Simulation](https://arxiv.org/abs/2007.14720), [DR2 BGS paper](https://arxiv.org/html/2507.01593v2), [EDR release](https://data.desi.lbl.gov/doc/releases/edr/vac/uchuu/). |
| Rosella / P-Millennium | 542 Mpc/h; 5040^3; particle mass about 1.06e8 Msun/h | One simulation volume | BGS population at z~.2; SHAM luminosities, formation-history colours, halo/subhalo information. Data on request; no locally verified full survey/field package | Very relevant independent galaxy–environment benchmark and possible pretraining source. Colour assignment is modelled physics, not observational truth. [Paper](https://academic.oup.com/mnras/article/505/1/325/6270909). |
| GLAM BGS | Published BGS demonstration uses 500 Mpc/h boxes; exact resolution of local v2 not certified | GLAM supports ensembles, but the number of independent phases underlying this BGS release was not established | BGS HOD from Uchuu described publicly. Local tested product truncates above z=.50; matched matter-field access unverified | Retain as a selection/clustering control; do not borrow the 2,000-mock count from a different ELG release. [BGS construction presentation](https://indico.in2p3.fr/event/29646/contributions/130834/attachments/80869/118930/Mocks%20for%20the%20DESI%20BGS%20.pdf). |
| MillenniumTNG | 500 Mpc/h hydro and matched 4320^3 DMO; larger neutrino simulations also exist | Flagship DMO pair uses fixed/inverted phases; hydro is one realization | Hydro/SAM populations and lightcones are scientifically attractive. Verified ready BGS+matter download not found; official site still advertises a planned 2024 release | Promising physics/transfer comparison, not an immediately available independent-phase BGS replacement. [Project](https://www.mtng-project.org/), [simulation paper](https://academic.oup.com/mnras/article/524/2/2556/7226466). |
| FLAMINGO | High-resolution hydro L1_m8: 1 Gpc, baryon mass 1.34e8 Msun, DM mass 7.06e8 Msun; flagship 2.8 Gpc has coarser mass resolution | Feedback/cosmology variations and observers are not an independent-phase ensemble; an inverted-phase DMO counterpart exists | Public data release; retained particle snapshots include z=0,.1,.2,.3,.4,.5. No ready complete BGS observational pipeline verified | Useful external hydro/feedback test. Does not presently meet all ensemble+BGS requirements. Lengths/masses here are physical Mpc/Msun conventions, unlike the h-scaled rows above. [Specifications](https://flamingo.strw.leidenuniv.nl/simulations.html), [release](https://flamingo.strw.leidenuniv.nl/data.html). |
| Existing TNG300 | 205 Mpc/h, with archived dense galaxy graph and native field already available | One primary volume; subvolumes/resolution levels are not independent phases | All ingredients for cheap controlled density and target-scale experiments; no ready DESI selection pipeline | Best immediate mechanism test and pretraining pilot. Limited volume and galaxy-formation assumptions preclude using it as sole survey validation. [Public access](https://www.tng-project.org/data/). |

The local-product constraints above are inherited from the actual September
audits, not inferred from simulation marketing or the general existence of
BGS mocks. See [mock choice and tested support](p12a_mock_choice_20260928.md)
and [alignment execution results](p12a_alignment_execution_results_20260926.md).
They do not mean all possible future Uchuu/GLAM products have the same cutoff.
HOLI reached the tested high-z interval but used Loa-calibrated n(z) and
placeholder magnitudes; that is not independent galaxy–matter validation.

The apparent Abacus complexity is concentrated at the interface between the
gravity simulation, halo cleaning, independently populated galaxy branches,
periodic replications, survey coordinates and observation selection. Here
HALO_INDEX and BOX_INDEX refer to different filtered objects. Stable parent
IDs and a retained transform would make this much simpler. Changing gravity
simulations does not automatically preserve those records, and a BGS catalogue
that lacks them recreates the same problem elsewhere.

## Pretraining and the remaining discriminating experiments

High-resolution pretraining followed by Abacus adaptation is a sensible
hypothesis to test. Dense data can teach geometric features; subsequent
training on varied sparse observations can teach how uncertainty changes with
selection. Pretraining cannot reconstruct information absent from a sparse
observation, and it may import an overly confident dense-population prior.
Its benefit must be judged on held-out simulation truth and posterior coverage,
not by maximizing observed property gradients.

| Experiment | Concrete comparison | What it would isolate | Required acceptance evidence |
|---|---|---|---|
| Fixed target, varying tracer density | On TNG, rebuild the same encoder inputs at full, .3 and .1 density; keep query galaxies and truth fixed | Whether the information loss shown above also affects the trained predictor | Eigenvalue errors, class confusion and uncertainty by density; spatially separated train/validation regions with context buffers |
| Fixed tracers, varying truth scale | Derive R2 and R7 from the same matter field using the same deposition/normalization; fit identical models | Target smoothing versus tracer sparsity | Same-query comparison and resolution/convergence checks; threshold .2 held fixed unless explicitly redefined |
| Frozen GAT versus P12-A transfer | One already exposed Abacus phase, identical observed tracers/core rows, both independently verified truths | Method transfer and scale mismatch | 2-model × 2-truth-scale matrix, density-only baseline, proper scores/calibration; not an architecture-only claim |
| A: Abacus from scratch | Current architecture, same ordered-eigenvalue target and posterior head | Reference performance | Fixed training phases, observation draws, optimizer steps and validation protocol |
| B: dense pretraining → Abacus | Same architecture pretrained on high-resolution periodic data, then adapted on the same Abacus set | Benefit of dense geometric pretraining | Improvement over A on unseen phases without coverage degradation |
| C: density-augmented pretraining → Abacus | As B, but pretraining includes thinning, RSD and appropriate masks/selection conditioning | Whether domain randomization is necessary | Compare B/C at matched total training budget; test density and redshift strata |
| D: longer Abacus-only training | Match the total updates/compute of B/C | Extra optimization versus transferable information | B/C must improve beyond this control, not just A |
| External simulation transfer | Evaluate frozen variants on a verified Uchuu or Rosella mapping | Dependence on HOD/SHAM/hydro assumptions | Identical truth definition and selection; clearly separate independent physics from independent phases |

Use the **same downstream architecture** for A–D. The old GAT weights cannot
be loaded directly into a voxel U-Net; that would require a separate
distillation/representation study. Prefer ordered continuous eigenvalue
targets with a shared posterior contract; derive classes afterwards. Keep
mass, colour and SFR out of inputs for the primary property-validation test.
If they are later used as features, their own trends cease to be independent
validation outcomes.

Use three optimizer seeds initially, preserve the repository's existing
reserved-phase policy, and tune only on already exposed validation data.
Report paired errors, proper scoring rules, interval coverage, class
reliability and worst density/redshift strata. Recalibrate after adaptation;
do not inherit dense-domain posterior calibration. Whole phases, rather than
neighbouring nodes in one graph, must supply the final independent test.

These are proposed experiments, not completed training runs. The present
deliverable closes the diagnostic/research report while explicitly leaving
matched-model causal attribution and pretraining efficacy open. No production
model, catalogue or science redshift cut has been changed.
