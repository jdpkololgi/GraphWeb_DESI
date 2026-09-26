# Mock-to-LOA Alignment Investigation

**Status:** active investigation plan, 2026-09-26. GraphWeb_DESI owns observation
and catalogue work; Illustris owns simulation truth, training and the science
log. This document is the canonical population-alignment subplan of the
[P12-A VAC roadmap](../../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md).
It supersedes the execution ordering in the earlier
[population-test protocol](p12a_population_test_protocol_20260926.md), retaining
that document and all prior audits as evidence. New experiments remain planned
unless explicitly marked completed below. No production catalogue is modified
by adopting this plan.

**Collaboration scope:** the user is a full DESI member and this work is under
the DESI collaboration. Internal releases, catalogues, production documentation
and mocks are in scope as primary sources. Public-release status is descriptive,
not an access or scientific-eligibility gate. Search authorized internal products
first; distinguish a product not yet located from one unavailable to the user.

## 1. Objective and scientific contract

Determine why the mocks fail to reproduce the Loa BGS-BRIGHT population, quantify
which mechanisms explain the discrepancy, and qualify a revised forward model
for environmental posterior inference. Matching n(z) is a necessary diagnostic;
qualification also concerns clustering, galaxy properties, observation losses
and the relationship between galaxies and matter.

The primary target stays **Loa/DA2, LSS v2.1, full BGS_BRIGHT, 0.15<=z<0.55**,
with the current ZWARN=0, DELTACHI2>=25 and SPECTYPE=GALAXY baseline. Keep a
separate DELTACHI2>40+GALAXY control. Preserve the actual PHOTSYS targeting
branches, including northern magnitude limits and historical recovery rules.
BAO absolute-magnitude-selected samples, ANY, intrinsic parents, observed full
catalogues and weighted clustering catalogues are separate comparison populations.
A sample/range change must be a named alternative, never an unnoticed fix.

The existing 5,436,413-row full-survey VAC is **provisional**. Production
completion is not science-release qualification. P12-A targets joint ordered
per-galaxy eigenvalues at its pinned epoch, smoothing and coordinate conventions;
alternative multi-epoch simulations do not automatically have the same estimand.

Two outputs are required:

1. An explanatory account: which upstream/observation differences are established,
   their paired effects on the residual, interactions and what remains unknown.
2. A validated practical solution: either an accepted revised mock/observation
   model with qualified posteriors, or a clear finding that the tested candidates
   do not support the requested VAC domain. An empirically calibrated selection
   may be useful without proving the unique physical cause; label that distinction.

## 2. Starting evidence: do not restart resolved checks

| Established result | Implication |
|---|---|
| Common-sky high-shell Loa/original-mock ratios1.856 SGC,1.839 NGC | Large discrepancy needing explanation; it grows more sharply in the finest high-z bins. |
| All256 conjunctions of eight cuts leave1.832/1.718 | Those extra flags cannot account for it; do not repeat without changed inputs. |
| Full-BRIGHT Loa v1.1/v2.1 broad counts differ<0.5% | This catalogue-version switch does not solve the mismatch. |
| Historic preparation replays after restoring Y5/IDs; tile-only change tiny on common sky | IN_Y flags alone are not the explanation. |
| Broader Loa bright tail at fixed z; passband equivalence unresolved | Luminosity support/evolution and observer photometry remain plausible. |
| Initial49 magnitude-only trials improved broad counts but failed fine-z/colour checks | No accepted repair; use as prior sensitivity evidence, not a production model. |
| ANY-02 matches aggregate BAO counts below0.4 but overpredicts at0.4–0.5; LoaBRIGHT-02 remains deficient | BAO workaround is evidence and a method to test, not a universal replacement. |
| Raw Abacusv1 ph000 and public Uchuu SV3 BGS are available | Test these alternatives before dismissing them. Their parity is not established. |

Sources: [joint cuts](p12a_joint_flags_magnitude_20260926.md),
[test phase](p12a_mock_test_phase_20260926.md),
[reference-mock review](p12a_reference_mock_review_20260926.md),
[selection closure](p12a_selection_closure_20260926.md), and their compact evidence.

## 3. Data releases and mock inventory

Live bounded inventory checked2026-09-26; see
[evidence](evidence/p12a_alignment_investigation_20260926/RELEASE_INVENTORY.json).
Directory existence does not establish production completeness or qualification.
Internal availability is sufficient for this collaboration investigation; the
public inventory is supplementary and must not constrain candidate selection.

| Data/product | Verified availability | Associated mock position and use |
|---|---|---|
| Loa/DA2 | Current target, internal full catalogue and LSS versions | Existing SecondGen Abacus BRIGHT/ANY products; primary development comparison. |
| DR1/iron | Public release and local Y1/iron LSS | Public Abacus/EZmock products documented; use matched sample and common-area/TARGETID checks, not an assumed independent survey. |
| EDR/SV3/fuji | Public release and local SV3 LSS | Readable Uchuu BGS ensemble for its SV3 geometry; useful first external population check. |
| DA3/matterhorn-v2 | Internal `DA3/LSS/matterhorn-v2/LSScats/v0/BGS_BRIGHT_full_HPmapcut.dat.fits`,11,686,801 header rows, including non-science/unobserved rows | No `DA3/mocks` directory in the checked root. Associated BGS mock package not identified; search production documentation/explicit linked locations before claiming absence. |
| Public DESI DR3/DR4 | Not listed in current official release index; corresponding local public directories absent | No public DR3/DR4 BGS mock release verified. Internal DA3 work is not proof of a public DR3 release. |
| DA4 | Absent from inspected survey/catalogs root | No associated BGS mock product identified in this bounded search. |

Official sources: [release index](https://data.desi.lbl.gov/doc/releases/),
[DR1 mock documentation](https://data.desi.lbl.gov/doc/releases/dr1/),
[Uchuu VAC](https://data.desi.lbl.gov/doc/releases/edr/vac/uchuu/),
[LSS matterhorn/DA3 configuration](https://github.com/desihub/LSS/blob/main/py/LSS/globals.py).
The local public/dr2 README says it is a placeholder: do not infer a public
spectroscopic release from the directory name. Distinguish DESI releases from
Legacy Imaging, Gaia or GAMA releases encountered in searches.

Candidate mock families:
- Raw Abacusv0/v0.1/v1 at `cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/`.
  v1 currently has only ph000;74.75m total rows versus63.90m in v0.1, with a
  different halo/coordinate schema and no saved IN_Y flags. Audit before use.
- Current Abacus BGS_ANY raw parent and processed '-02' controls. Preserve
  original magnitudes/IDs and distinguish Kibo versus Loa selection.
- Public Uchuu `public/edr/vac/edr/uchuu/v1.0/BGS-BRIGHT_Uchuu/`.
 102 small-footprint realizations are not102 independent full-Loa boxes. Delivered
 columns contain photometry/K/E but no halo IDs. Locate the actual DR2 BGS parent
 and particle/halo linkage separately; previously searched DA2 locations did not
 provide a readable full-Loa BGS product.
- GLAM/Holi/other families: directories or LRG catalogues alone do not establish
 a BGS alternative. Inventory tracer, recipe, truth access and resolution before
 selecting a candidate. Covariance mocks need not support environmental truth.

## 4. Work packages, dependencies and decision gates

### W0 — Freeze comparisons, provenance and exposure ledger (first)

Write a machine-readable manifest for each actual input: simulation/cosmology,
phase, generator revision/parameters, magnitude bands and zero points, K/E
conventions, h units, footprint, targeting version, assignment/masks, redshift
quality, weights, schema, row IDs and file fingerprint. Separate reconstructed
recipes from executed commands supported by production logs/headers.

Record already-exposed ph000 and ph002–006 roles from the authoritative ledger.
Do not open reserved phases, even for a cheap distribution check, until their
use is authorized by that ledger and the frozen confirmation contract. Existing
exposures cannot become independent holdouts by relabelling.

Freeze common angular support, cap/PHOTSYS subdivisions, dz=.01 diagnostics and
four primary shells. For new parents regenerate footprint membership from the
same tile/mask definitions. Use randoms to assess holes/subpixel completeness;
occupied-pixel matching alone is an initial control, not exact angular parity.
Keep counts per area separate from n(z) per volume and document cosmology/units.

**Deliverables:** input manifest, sample crosswalk, phase ledger, baseline plots.
**Gate G0:** meaningful comparison contract exists; unknown observables remain
explicitly unknown rather than assigned invented mock values.

### W1 — Fast alternative/version screens (after G0; independent of LF fitting)

1. Compare Abacusv1 andv0.1 on the same exposed ph000, same geometry and physical
   apparent-magnitude definition. Record parent support, n(z), r/colour versus z,
   bright-tail bounds, central/satellite mix and halo-ID/coordinate compatibility.
   Determine what changed in the generator; more rows alone is not a pass.
2. Compare public Uchuu BGS with its intended SV3 data, identical cuts and
   completeness treatment. Use all eligible realizations for descriptive scatter,
   accounting for shared simulation volume. Compare intrinsic-to-corrected or
   observed-to-observed populations; never intrinsic counts to raw Loa successes.
3. For a verified full-sky/DR2 Uchuu parent, construct the Loa footprint and
   observation comparison. Passing SV3 does not establish passing Loa.
4. Prioritize internal DA3/matterhorn and its production/mock provenance alongside
   DR1 controls. Inspect internal LSS configurations, production logs and explicitly
   referenced mock locations; do not stop at a missing DA3/mocks directory. Record
   exact version, observation interval and target/mask definitions. Compare DR1
   and DA3 with Loa on common TARGETIDs and sky, then separately the
   added area/epochs. Separate changes in redshifts, photometry, quality and
   footprint. Check internal version/blinding status before interpreting positions.
   A later catalogue is a diagnostic, not a silent replacement of the science target.

**Deliverable:** candidate scorecard with explicit support/truth/provenance limits.
**Gate G1:** retain promising candidates; reject unsupported comparisons, not
whole simulation families. If a newer Abacus recipe fixes parity, reproduce that
mechanism before investing in bespoke tuning.

### W2 — Reconstruct the baseline observation chain and close remaining cuts

Replay raw parent -> targeting -> masks -> assignment -> spectral success ->
science sample. Match historical baseline on exposed phase before altering it.
Retain the previously passed flag audits; finish unresolved photometry-at-targeting,
mask/exposure/bad-fibre versions and accepted fibre-status semantics.

Pin the absolute-to-apparent mapping, pivot redshift, distance/h convention,
colour-dependent K-correction, existing luminosity remapping and stored E terms.
Resolve whether a magnitude already contains evolution before applying another
term. Compare same-galaxy Legacy/SDSS photometry as a function of colour, z and
PHOTSYS, using clean joins and explicit selection. Use appropriate FastSpecFit
K-corrections and completeness-aware LF estimation; do not fit an LF to raw
success counts without the selection function.

Mocks lacking DELTACHI2/SPECTYPE/fibre magnitude need a justified observation
model, not numerical cuts on nonexistent quantities. True z of real failures is
unavailable: constrain success with repeat/deeper spectra or observable-based
models, carry uncertainty, and separate this from assignment completeness.

**Deliverables:** replay report, mapping residuals, stage retention curves,
resolved/open provenance table. **Gate G2:** exact/quantified baseline agreement
and no unresolved convention capable of masquerading as the proposed effect.

### W3 — Controlled mechanism tests (after G2)

Use identical halos, geometry and seeds. Start with one exposed phase for speed;
replicate promising effects on other exposed phases. Test:

| Contrast | Fixed quantities | Question |
|---|---|---|
| P-only changes | Q, K and photometry | Is density evolution sufficient? |
| Q-only changes | P, K and photometry | Is luminosity evolution sufficient? |
| K/passband changes | P/Q and intrinsic population | Does observer-band mapping explain the redshift/colour dependence? |
| Bright-tail/LF-shape changes | Explicit fixed evolution model | Is the parent luminosity support or shape wrong? |
| Joint P/Q/K model | Same geometry and observation recipe | Are interactions required? |
| Observation-success variants | Fixed parent photometry | How much can selection/failure modelling contribute? |

Baseline nominal P1.8/Q0.7 must be verified against the executed generator.
Vary its actual LF/remapping, not a final plotting E-correction or arbitrary
weights. Begin with small symmetric perturbations to measure responses; freeze
physically justified ranges from LF constraints before a wider fit. Document
when deeper parent regeneration is required by finite magnitude/halo support.

Measure effect sizes and residuals in each shell/cap, conditional r/colour
quantiles and luminosity functions. Quantify interactions rather than summing
non-independent effects. A successful ablation establishes sensitivity; a causal
explanation also needs baseline provenance and correct joint predictions.

**Deliverable:** paired response matrix and ranked mechanisms, with rejected
hypotheses retained. **Gate G3:** candidate explains counts without worsening
other required observables beyond the registered tolerances.

### W4 — Empirical ANY selection candidate (alongside W3 after G2)

Test the user's proposed smooth redshift-dependent absolute-magnitude cut on
the broad parent, without assuming the BAO coefficients apply to full BRIGHT.
Verify there are enough eligible parent objects in every required shell; no cut
can add objects that do not exist. Distinguish raw parent depth from successfully
observed FAINT targets.

Fit M_lim(z), or a smooth probabilistic selection if needed, against Loa's
expected selection density, with explicit absolute-M conventions and estimated
observation losses. Regularize the function and constrain model complexity.
Do not reproduce every observed radial fluctuation: that would absorb real
structure. Reserve NGC as a transfer check for an SGC-tuned model, with PHOTSYS
calibration differences explicit; both caps remain development data.

Keep original apparent magnitudes and added selection fields distinct. If the
candidate admits numerically faint objects, label it an empirical tracer mapping;
it does not establish physical passband equivalence. Compare it fairly with the
LF-based candidates using the same diagnostics and observation pipeline.

**Deliverable:** reproducible candidate selection and support/complexity report.
**Gate G4:** sufficient support and promising joint population predictions; count
agreement alone is not permission to replace inputs.

### W5 — Fresh observation processing of shortlisted candidates

Run current validated DESI preparation/target priorities, masks and assignment
for changed target populations. Include competing populations relevant to fibre
assignment. Do not inherit FAINT assignments or old BRIGHT retention fractions.
Use measured/validated spectral-success treatment and bracket uncertainty where
mock observables are absent. Retention approximations remain screening-only.

Rebuild randoms, response/expected-density fields and survey feature inputs for
the chosen sample. Check pinned coordinate conventions, observational RSD versus
physical offsets, and source/selection consistency. PIP/FKP weights are used
only for the estimators they support; they do not create missing galaxy positions.

**Gate G5:** reproducible observed-catalogue replay, stage-loss closure and no
input/model convention mismatch. Archive commands, seeds, source hashes and logs.

### W6 — Joint statistical validation and candidate freeze

Evaluate on data not used for each tuning step, preserving spatial correlations.
Compare raw observed counts, weighted population estimates and intrinsic-parent
quantities separately. Mandatory diagnostics:

- Broad and fine n(z), sky/PHOTSYS dependence, r–z and colour–z distributions,
  conditional tails and fibre magnitudes where physically defined.
- Luminosity-dependent w_p, redshift-space monopole/quadrupole, and neighbour or
  counts-in-cells statistics on scales relevant to the encoder and7Mpc/h target.
- Simulated host mass, central/satellite fractions and tracer–matter relation;
  consistency with justified external galaxy–halo constraints, not assumed truth
  for real galaxies. Compare alternate priors where these remain weakly known.
- Completeness/response closure and high-z feature-domain overlap.

Use phase scatter, spatial resampling and observational uncertainties; do not use
independent-galaxy Poisson errors as the complete uncertainty. Five exposed phases
cannot support inversion of a high-dimensional sample covariance. Reduce the
predeclared summary dimension or validate a suitable covariance model. Report
practical effect sizes as well as tests, with multiple comparisons accounted for.

**Provisional screening tolerances:**5% broad-shell and10% dz=.01 count residuals
where uncertainty supports that precision; investigate coherent excesses even
when totals match. These are not release criteria and must not be relaxed after
seeing an unfavorable candidate. W0/W2 must freeze final joint tolerances using
measurement uncertainty and downstream sensitivity before selection/confirmation.
When data cannot distinguish models, propagate that uncertainty rather than
calling a candidate validated. Fitted summaries are not independent validation.

**Gate G6:** freeze source, parameters, sample, diagnostics, tolerances and
confirmation exposure. A candidate can pass population checks without uniquely
identifying the historical cause; report both conclusions separately.

### W7 — Environmental inference and VAC replacement gate

First assess frozen P12-A feature/summary residuals and conditional mock coverage
under the candidate; this measures transfer and the need for retraining.
If the revised population changes the learned summaries or conditional posterior,
regenerate OOF summaries, refit the posterior and retrain encoder weights where
required by evidence. Do not repair a changed forward model with calibration
alone. Preserve truth epoch, smoothing, halo/particle mapping and ordered targets;
Uchuu requires a fresh truth-compatibility audit.

Use the existing registered coverage/TARP and conditional-calibration protocol,
sharpness/bias checks, and selection-sensitivity tests. Perform confirmation only
on eligible reserved simulations after freezing decisions. Real-data closure
cannot prove real-DESI frequentist coverage without truth; document simulation
and observation-model uncertainty explicitly.

**Gate G7:** independent mock confirmation, golden replay, bounded Loa inference,
input-distribution acceptance, then restartable full-survey replacement production.
Version the replacement and compare it with the provisional VAC. Science release
requires these gates; a matched n(z) or a completed Slurm job is insufficient.

## 5. Execution order, outputs and stopping rules

Immediate sequence: W0 -> W1 screens and W2 replay -> W3/W4 candidates -> W5 ->
W6 -> W7. Independent preparation may proceed concurrently without opening
confirmation data. A failed prerequisite stops only dependent work.

| Stage | Current status | Next concrete artifact |
|---|---|---|
| Prior cuts/version/BAO audits | Completed, limitations retained | Existing reports/evidence linked above |
| Later-release and alternative inventory | Header/directory checks completed | RELEASE_INVENTORY.json and earlier ALTERNATIVE_INVENTORY.json |
| W0 | To freeze | Manifest, exposure ledger and registered comparison contract |
| W1 | Not run | Matched Abacusv1/ph000 and Uchuu/SV3 population scorecards |
| W2 | Partially closed by existing audits | Exact executed generator/mapping replay and remaining observation closure |
| W3/W4 | Not run under this protocol | Factorial mechanism results and full-BRIGHT ANY candidate |
| W5–W7 | Dependent, not launched | Observed candidates, qualification reports and versioned VAC decision |

Keep source/configuration/scientific decisions and compact evidence in Git; large
catalogues/checkpoints on Scratch; CFS inputs read-only. Use cosmic_env (or its
approved mirror), graphify from that environment, compute nodes for substantial
scans, the two-allocation limit and existing compute authorization. Log failures
as well as successes. Freeze production source and use restartable bounded jobs;
no automatic unlimited sweeps or repeated retuning of confirmation data.

Stop and report a decision if: all parents lack required support; the required
production provenance/halo truth cannot be recovered; selection changes need a
new scientific sample; all candidates fail joint validation; or remaining
uncertainty prevents release. Options include a new generator/deeper parent,
a separately qualified restricted-domain VAC, or continued provisional status.
Do not switch to DR1/DA3, lower z_max, or choose a simulation suite merely because
it gives a more convenient count curve.
