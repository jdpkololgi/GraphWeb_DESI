# GraphWeb field-posterior inference: fresh-eyes briefing

Status: 16 October 2026. This is a methodological handoff, not a claim that the
field models are calibrated. It distinguishes the simulation evaluation from
the observed-DESI visualization programme. All numbers below refer to the
specified completed evaluations, not a new experiment.

## 1. Scientific aim and the two different posterior products

We want to infer the cosmic matter distribution and gravitational tidal
environment from the sparse, survey-selected positions of observed galaxies.
Galaxies are biased, incomplete tracers of matter; survey boundaries, radial
selection and unobserved large-scale matter make this an uncertainty problem,
not just a deterministic image-reconstruction problem.

The tidal tensor is the Hessian of the gravitational potential. Its ordered
eigenvalues are lambda1 <= lambda2 <= lambda3. Counting eigenvalues above the
threshold 0.2 gives void/wall/filament/cluster classes (0/1/2/3). Those classes
are derived summaries, not the primary learned field target. This is T-web,
not velocity-shear V-web.

There are two distinct products:

- **P12-A: per-galaxy eigenvalue posterior estimator.** A coordinate-aligned
  conditional flow estimates the three ordered eigenvalues at each galaxy.
  Its registered blind mock evaluation passed its ten release gates on
  4,897,905 galaxies. This is conditional-on-mock/HOD evidence, not DESI
  validation, HOD marginalization, exact conditional calibration, or a spatially
  joint field posterior. Independently sampling galaxies does not make one.
- **E2E: field-level generative inference.** Sample a positive matter-density
  field, then apply a fixed physical tidal operator and eigensolver. A draw
  supplies mutually related density, tensors and environments within its
  modeled domain. Current finite-domain or independently tiled products are
  not a qualified jointly conditioned full-survey posterior.

The current E2E fits are voxel-based conditional generative models, not the
historical galaxy-graph GAT classifier or graph-NPE pipeline. “End-to-end” here
does not mean the entire coarse/fine hierarchy is trained jointly: its stages
are trained separately, and the coarse model is frozen for the newer fits.

## 2. How many posterior draws were actually used?

| Product or evaluation | Actual ensemble | Interpretation |
| --- | --- | --- |
| Earlier mixed-checkpoint CFM Abacus panel (September 26) | 32 coarse–fine draws per region pair and seed; 64 pairs x 2 seeds = 4,096 primary draws | Not four; extra 8-draw subsets assessed NFE refinement |
| October 6 legacy-CFM Tier-1 development analysis | 32 draws per registered pair, 32 pairs / 64 cores in ph012/013 | Regional and marginal findings use a 32-draw reference |
| Latest expanded U-Net / transformer Abacus evaluation (October 14 closeout) | 32 draws per model/seed/region/coarse-mode case, NFE128 | 768 primary cases across both architectures, seeds, six phases and two coarse modes |
| Latest epoch-24 versus epoch-32 comparison | 32 draws per case; 64 paired cases per architecture | Geometry-selected development subset, not new confirmation |
| Observed Loa two-tile display pilot (October 15) | Four draws per fit/tile; six fits x two tiles = 48 fine draws | Exploratory display; observed data have no matter truth for calibration |
| Loa engineering preflight (October 16) | Four draws across 16 tiles, with two tiles extended to 32 | Throughput/replay/physical checks, not calibration |
| Newly approved survey run | Four-draw preview, then 32 only after review | Real-data inference, not a replacement Abacus calibration test |

Crucial distinction: the latest **sampled-coarse** evaluation uses 32 different
coarse draws, each with one addressed fine-noise draw. It samples the composed
law q(C|G) q(U|C,G); it is not 32 fine samples conditional on a single fixed C,
nor 32 x 32 nested draws. Shared reference files store four draws per chunk,
and the reader selects chunk `4*floor(draw/4)` and within-chunk index `draw%4`.
Eight such chunks supply 32 coarse draws. Four is a storage batch size there.

A direct metadata recount on16October found900 saved evaluation bindings, all
with32 draws:448 sampled and448 oracle at NFE128 (primary plus checkpoint
stability), and two of each at NFE256. The896 NFE128 addresses comprise768
primary final-checkpoint cases and128 stability cases. This recount checks
the saved ensemble contract; the existing closeout supplies the completed
payload audit. It is not a new read of every large field array.

The **true-coarse/oracle** mode holds C at the simulated truth and generates
32 independently addressed fine-noise draws. It tests a different conditioning
problem. The sampled mode marginalizes approximate coarse uncertainty; the
oracle mode does not. One fine sample per independently sampled C is a valid
Monte Carlo construction for the composed marginal law; it does not separately
resolve the conditional distribution at each fixed generated C.

Two training seeds are two fitted estimators, not twice as many independent
simulation universes. Matched noise between architectures/checkpoints makes
comparisons paired; it is not extra independent evidence. NFE128 means 128
neural evaluations per sampled trajectory, not 128 posterior samples.

## 3. Simulation and observation processing

Training uses AbacusSummit base c000, a fixed cosmology and HOD/mock selection,
with a 2,000 Mpc/h periodic simulation box and matter snapshot at z=0.2.
The galaxy observations are BGS_BRIGHT cutsky/alt-MTL mocks with survey response.
The science coordinates span observer redshift 0.15–0.55, with observation
context 0.10–0.60 and the recorded sentinel exclusion. This does **not** make
the matter training target an evolving lightcone or marginalize cosmology/HOD.

Matter processing uses field and halo particle subsamples A+B (3%+7%), TSC
assignment to a 2048^3 native lattice, full-box mean normalization and one
Gaussian smoothing with R=7 Mpc/h. Native spectral tidal references are kept
independently of the learned field. Regional mean density is preserved; regions
are not independently renormalized to unit mean. Coordinate conversions use
h=0.6766 and recorded grid origins/alignment. Density and tensor interpolation,
finite-domain closure and discretization are separate error sources.

Observed galaxies enter through counts and survey-response features, not matter
truth, target eigenvalues, phase IDs or filenames. Successful-redshift selection
requires finite positive Z_not4clus and ZWARN=0. Angular/random response and
radial selection are explicit conditioning information. The observation base
lattice is 5 Mpc = 3.383 Mpc/h; fine fields use 6.766 Mpc/h cells and wide coarse
fields use 27.064 Mpc/h cells. Inherited selection transforms are recorded;
new model normalization is fitted only on the thirteen training phases.

The twelve local channels are counts, random support, angular response,
apodized exposure, expected counts, log count ratio, distance to support
boundary, sampling intensity ntilde, three line-of-sight components, and
observer redshift. The twelve wide channels are counts, expected counts,
random support, angular response, apodized exposure, ntilde, geometry-valid
fraction, log count ratio, three line-of-sight components, and observer radius.
These features encode selection; they do not constitute a separately calibrated
explicit likelihood or prove transfer to real DESI.

The original fine domain has 64 x 48 x 48 cells. The expanded models use
96 x 64 x 64, with the same 48^3 wide coarse context (~1,299 Mpc/h per side).
Both score the same two adjacent 16^3 inner cores (8,192 voxels altogether).
Geometry selects cap, redshift-shell and boundary/interior strata without
choosing targets or favorable predictions. Contexts may overlap even when
owned cores do not; regions and voxels are not independent observations.

The thirteen training phases are ph000,002,003,007–011,020–024. Expanded fits
use 48 domains per phase, 624 total. Their voxel exposure per epoch equals
the older 1,664 smaller domains. ph012/013 are development; ph014–017 were
excluded from these fits but exposed in earlier analyses, so are out-of-fit
replication, not fresh blind confirmation. ph018/019 remain untouched.
P12-A's ph001 blind opening is not permission for E2E to use it.

## 4. Probabilistic factorization and physical constraints

Let G be observation/selection channels, rho the positive R7 matter density
in units of the cosmic mean, C its 4^3-cell block means, and U the block-centred
log-density shape. The modeled factorization is

    q(C,U | G) = q_coarse(C | G) q_fine(U | C,G).

The coarse variable is standardized log C. The fine variable is standardized
block-centred log rho. Projection subtracts each block mean. The decoder is

    rho_i = C_block * exp(U_i) / mean_block(exp(U)).

It uses a stable block softmax, so rho is positive and every block has exactly
the supplied coarse mass up to numerical precision. These are structural
identities, not learned calibration evidence. Coarse uncertainty cannot be
repaired by redistributing fine density inside an exactly conserved block.

Fine training supplies the paired **true** coarse field. Real inference supplies
a generated coarse field. Errors/shift in this interface can matter, but an
oracle improvement alone is not an additive decomposition of error, and does
not prove that the interface is the sole cause. Pairing a training residual
with arbitrary independent coarse draws would change the target joint law.

## 5. Objectives, architectures and numerical sampling

The current six-fit Loa comparison is three fine configurations x seeds17/29:
the earlier joint CFM velocity model; expanded U-Net clean-endpoint CFM; and
expanded transformer velocity CFM. Each seed's frozen coarse CFM is shared
across the three fine models. Legacy checkpoints are coarse step13,312 and
fine step26,624; expanded fine checkpoints are step9,984 (32 epochs), using EMA.

For conditional flow matching, Gaussian noise epsilon and target x define the
linear bridge z_t=(1-t)epsilon+t*x and target velocity x-epsilon. Fine noise,
state and velocity live in the block-zero subspace. Minimize Fourier-weighted
velocity squared error, using a training-only power profile. Amplitude weighting
is proportional to P(k)^(-0.30), with recorded normalization/floor; because it
is squared in the loss, the squared-error weight is P(k)^(-0.60). This is a
deliberately weighted finite-capacity regression objective, not a direct
coverage loss. Low bridge loss does not establish posterior calibration.

The expanded U-Net has base width24, three resolution levels, time conditioning,
skip/residual paths, bottleneck self-attention and cross-attention to local/wide
conditioning with physical positions (~8.40M parameters). It predicts a clean
endpoint, converted to velocity by (xhat-z)/(1-t), with denominator clamped at
1/1024. Endpoint conversion makes near-final integration behavior important.

The expanded transformer has eight width240, six-head blocks, 4^3-cell patches,
physical position/type encodings and time modulation (~8.94M parameters).
Its 6,144 fine and 1,728 wide tokens interact; it predicts velocity directly.
Backbone and parameterization differ simultaneously, so this is not a clean
isolated test of “CNN versus transformer.”

Expanded training uses batch2, AdamW, cosine learning rate 1e-4 to 1e-5,
gradient clipping0.5 and EMA0.999. Forward computation is BF16; parameters,
optimizer state and loss are float32. Data/noise addresses are matched between
arms. These fits have 9,984 updates; equal exposure does not guarantee equal
optimization or convergence. Earlier diffusion/denoising ablations are historical
diagnostics, not the objectives of these latest three CFM configurations.

Sampling integrates the learned flow with Heun, 64 two-evaluation steps at
NFE128. The legacy sampler's argument64 denotes steps, also NFE128. Matched
NFE256 reruns change numerical resolution, not the number of independent draws.

Density is mapped through a fixed potential-Hessian operator. For a periodic
native reference, T_ij(k)=(k_i*k_j/k^2)*delta_R7(k), DC zero; smoothing has
already been applied once. Eigenvalues come from the resulting symmetric
tensor, not three unconstrained network outputs. Trace, ordering, positivity
and block mass are checked independently of statistical scores.

The latest Abacus primary evaluation retains a frozen reference full-cap
exterior, replaces the generated local density and applies the free-space
tensor difference plus the same conservative reconstruction correction.
Zero replacement must recover the reference. This is a modular replacement,
not a newly jointly sampled exterior conditional on the replacement.

Loa currently uses a different, explicitly labeled finite local-plus-wide
free-space closure with mean density beyond the wide domain. Tiles have
independent noise addresses. A cap-wide display of tiles is therefore not a
sample from one coherent full-cap posterior. Do not pool its scores with the
corrected Abacus full-cap-reference evaluation.

## 6. What is actually concerning about calibration?

The latest exposed Abacus evaluation contains four phases x16 region pairs
x2 training seeds =128 cases per architecture/coarse mode, but only64 distinct
region pairs in four simulation phases. Density CRPS is averaged over correlated
voxels; regional energy scores a six-vector: density, tidal spread
(lambda3-lambda1)/sqrt(2), and curvature (lambda1-2lambda2+lambda3)/sqrt(6)
in each core. Lower scores are useful, not sufficient for joint calibration.

| Latest sampled-coarse result | U-Net endpoint | Transformer velocity |
| --- | ---: | ---: |
| Density empirical CRPS | 0.134611 | 0.177831 |
| Regional empirical energy score | 0.063764 | 0.111060 |
| Interpolated 5%–95% voxel coverage (not exact 90% reference) | 0.875940 | 0.688046 |

The frozen reference regional energy is0.065790: U-Net improves by3.08%, not
uniformly across summaries/shells. Transformer regional tidal spread is about
21% low in both cores, with negative mean error in every one of the128 cases.
The bias persists with true coarse input and greatly exceeds the measured mean
closure-on-truth discrepancy. Small first-case NFE changes do not repair it.
This is evidence of a substantial structural problem, not merely a tiny
coverage deviation. It still does not identify a unique architectural cause.

U-Net has useful reconstruction and much smaller discrepancies, but regional
ranks remain uneven. From epoch24 to32, sampled-coarse density CRPS improves
0.69% while regional energy worsens1.17% on the development stability subset.
Some matched draws shift materially. Its first-case oracle NFE128/256 change
reaches0.929 posterior SD for one regional component; numerical convergence
has not been demonstrated panel-wide. “Not yet qualified” is justified;
“intrinsically incapable of calibrated inference” is not.

Earlier legacy-CFM results were also mixed, not universal failure: development
cellwise density/eigenvalue marginals passed the recorded tolerance screens,
while regional paired-gap coverage was56.25/68.75% against an87.879% finite-rank
reference. These findings concern particular summaries, operator and panel.
They do not mean every per-galaxy posterior or downstream relationship is wrong.

## 7. Monte Carlo limitations and defensible next checks

The concern about four draws is valid in principle, but it is not the ensemble
size of the principal Abacus results above. Even32 is modest:

- For iid draws and finite variance, posterior-mean Monte Carlo SE is
  sigma/sqrt(M): 0.50sigma at4, 0.177sigma at32, 0.088sigma at128.
- A class-frequency estimate has SE sqrt(p(1-p)/M), at most0.25 at4 and0.088
  at32. These are sampling errors conditional on a model, not guarantees that
  its probabilities are correct. Spatially correlated galaxies cannot simply
  be counted as independent trials to shrink an error bar.
- Under Gaussian iid sampling, relative sample-variance SE is sqrt(2/(M-1)):
  about82% at4 and25% at32. Non-Gaussian tails can be worse.
- With M continuous exchangeable posterior draws and truth, the truth rank is
  uniform over M+1 positions. The saved central event of zero-based ranks2…30 at M=32
  has probability29/33=87.879%, not90%. At4 even the sample min–max interval
  has unconditional reference coverage3/5=60%. Interpolated5%/95% quantiles
  have a distribution-dependent finite-M reference; do not mechanically
  compare them with90%, or substitute29/33 as their exact null.
- One truth per observation cannot establish exact conditional calibration.
  Region/phase dependence and repeated scientific inspection limit population
  claims even if Monte Carlo error becomes negligible.

The latest scorer uses empirical CRPS/energy V-statistics (pairwise denominator
M^2, including the zero diagonal). For estimating scores of the underlying
continuous predictive law, the fair/U-statistic version excludes self-pairs
and uses M(M-1). At equal M the difference is still spread-dependent, so small
score gains should be checked for estimator sensitivity. The closeout's
independent arithmetic verifies the implemented empirical scores; it does not
remove finite-ensemble bias. No fair-score rerun is claimed here.

A useful follow-up, **not launched by this briefing**, is a frozen geometry-
selected mock panel with independently addressed larger ensembles, reporting
4/8/16/32/64/128-prefix behavior and replicate-noise variability. Choose draw
counts from the needed precision, not128 as a magic threshold. For a separate
fine-at-fixed-coarse question, use nested repeated fine draws at selected C and
report that different estimand. Use discrete/tie-aware ranks, fair scores where
appropriate, paired numerical refinement, and phase-level uncertainty without
treating millions of voxels as independent universes. More draws cannot repair
a biased model, add independent phases, or validate transfer to DESI.

The approved Loa survey run is not this additional mock experiment. Its preview
is an engineering/coverage review, and its32-draw products remain exploratory.

## 8. Questions for an independent reviewer

1. Which discrepancies remain after separating finite-M error, score-estimator
   bias, numerical integration and operator/exterior approximation?
2. Is the desired law the marginal field given G or fine structure at a fixed C,
   and are the diagnostic ensembles adequate for that specific question?
3. How should true-coarse training versus sampled-coarse inference be tested
   without changing the intended joint law through arbitrary augmentation?
4. Are spectral weighting, endpoint conversion and optimization causing the
   observed structure/uncertainty tradeoffs? What matched control would decide?
5. What can be learned with the few available independent simulation phases,
   and what calibration claim would require fresh phases or HOD variation?
6. Which downstream galaxy/environment estimands are robust to the measured
   field defects, rather than assuming all succeed or all fail together?

## Source map (local repository records, checked 16 October)

The shared `TNG/Illustris/SCIENCE_LOG.md` controls current status over old plans.
Primary explanatory records in that repository:

- `configs/e2e_coupled_data_v1.json`: data, phases, target and geometry contract.
- `docs/e2e_cfm_mixed_20260926.md`: earlier32-draw panel specification.
- `docs/e2e_cfm_tier1_and_levers_20261006.md`: marginal/joint evidence and scope.
- `docs/e2e_wholefield_exploratory_pilot_20261013.md`: expanded model/training law.
- `docs/e2e_wholefield_pilot_closeout_20261014.md`: final results and MC caveats.
- `workflows/sbi/e2e_wholefield_pilot_evaluate.py`: actual32-draw loop, four-draw
  chunk addressing, oracle distinction, scoring definitions.
- `workflows/sbi/e2e_wholefield_pilot_models.py` and
  `e2e_wholefield_pilot_train.py`: exact architecture, bridge, objective, sampler.
- `workflows/sbi/e2e_vdm_context_models.py`: shared projection and positive decoder.
- `workflows/sbi/e2e_coupled_condition_products.py`: exact observation channels.

Saved final evidence is
`/pscratch/sd/d/dkololgi/abacus/e2e_field_v3/wholefield_pilot_20261013_v1/CLOSEOUT_FINAL_20261014.json`
(768 primary cases; posterior_qualified=false). GraphWeb_DESI's
`docs/loa_model_page_pilot_20261015.md` and
`docs/loa_survey_preflight_20261016.md` document the distinct observed-data runs.
