# P12-A: magnitude term, n(z), and regeneration assessment

2026-09-26. Diagnostic-only continuation; no mocks regenerated, n(z) imposed,
training performed or VAC replaced. User requests committing and logging work.

## What the existing magnitude term can and cannot establish

The raw ph006 numerical identity remains `m - M_stored - DM - K ~= -.8(z-.1)`.
The author's public generator applies an evolving luminosity function and uses
its evolved M directly in m=M+DM+K. Its published LF parameter table Q=.7 does
not identify the additional stored-column convention. Targeted searches of the
retained LSS mock tools and public Abacus generator did not locate the actual
-.8 transformation. The author's smith2024 public repository inventory supplies
paper data, not a recovered execution record for this local v0.1 product.
No reserved phase data were read. Therefore the earlier suggestion that this
term might already compensate for photometric differences remains a hypothesis.

There are at least two observationally degenerate interpretations of this
column identity:

1. Stored M is corrected to a reference epoch; converting back to M(z) requires
   the -.8 term. This is a definition, not an extra observer-band adjustment.
2. Apparent magnitude was subsequently brightened at fixed physical M, perhaps
   as an empirical mapping. Its impact is already present in the flux selection.

Algebra on the saved paired-SED sample tests the second interpretation only:

| z | SDSS r − DECam r | Existing term | Median paired sum |
|---|---:|---:|---:|
| .15–.25 |+.093|−.080|+.010|
| .25–.35 |+.106|−.160|−.057|
| .35–.45 |+.132|−.240|−.096|
| .45–.55 |+.246|−.320|−.076|

The term is large enough to cancel or over-cancel the synthetic filter
contribution. This makes adding another filter offset particularly unjustified.
These are conditional calculations, not evidence of why the term was applied.
Petrosian/model-flux effects are not included; medians do not add exactly.
An explicit producer definition or a joined before/after catalogue with pinned
physical absolute magnitudes is still required to distinguish the hypotheses.

## Where n(z) enters the checked BGS path

Source snapshots/hashes are in docs/evidence/p12a_nz_magnitude_20260926/.
Both the home preparation script and installed LSS script were inspected; these
are candidate source evidence, not certified historical execution receipts.

- `mask_secondgen(nz=1, foot='Y1')` constructs a STATUS bitmask. In the relevant
  second-generation branch, dark tracers apply it; bright BGS sets idx_main=idx.
  The presence of the name nz in this file therefore does not prove BGS was
  selected using a supplied n(z).
- BGS BRIGHT membership is assigned by R_MAG_APP < rbandcut (observed products
  previously verify19.5). FAINT is selected separately. No tabulated BGS n(z)
  argument is consumed in that branch.
- A generic --downsampling option defaults to n. Its BGS branch refers to a
  downsampling dictionary defined only for dark tracers in this inspected
  version. It is not an operational arbitrary BGS n(z) interface as written.
  We have not established that this option ran or that this source is identical
  to the producer's executed script; no production bug is inferred from it.
- The public Abacus generator obtains its flux-limited radial population through
  the evolving LF, rest-colour model, K-corrections, cosmology and flux limit.
  It resizes luminosities through cumulative-LF matching. No external BGS n(z)
  input was located in its inspected top-level generator.
- P12's fitted ntilde(z) is an expected-density input derived from the training
  population; changing it does not add missing tracers to the observed mock.

The user clarified that they remembered the capability generally, not a specific
entrypoint. A broader source search **confirmed an explicit BGS implementation**:

- LSS scripts/mock_tools/calibrate_nz_prep.py defines get_nz_bgs(nzfile),
  downsample_aux_bgs and calibrate_nz_bgs. It reads columns0/3 from the input
  table and selects with `ran < n_target(z)/n_parent(z)`.
- prepare_holi_bgs.py calls it with the actual Loa BGS_BRIGHT full-HPmapcut n(z)
  path, then selects STATUS==1 as its bright sample. This is a concrete example
  of the capability the user recalls, in a separate Holi pipeline.
- Ratios above1 keep all available parent objects; no new galaxies are created.
  This cannot close a deficit in the already bright-selected Abacus parent.
  Selecting from a larger BGS_ANY parent is possible, but changes the meaning
  of BRIGHT membership unless a physical luminosity/targeting model justifies it.
- The inspected Holi helper also assigns placeholder absolute magnitudes to
  sample categories. It is not interchangeable with our photometric/HOD-labelled
  Abacus input just because its n(z) matches. Its suitability for environmental
  inference would require separate validation.

Both source files are saved and hashed in DOWNSTREAM_NZ.json. No Holi catalogue
was opened, no n(z) calibration was executed, and no statement is made that this
later/current code was part of our historical Abacus production.

## Could an absent or wrong n(z) cause the discrepancy?

Yes, an inappropriate radial-population/selection model can cause a discrepancy.
But distinct operations must be separated:

| Operation | Effect | Interpretation here |
|---|---|---|
| Calibrate LF/HOD/population versus z | Changes which galaxies and host populations exist | Could remedy an established parent-population failure; requires clustering and environment validation |
| Thin a sufficiently dense parent with p(z)<=1 | Removes galaxies | Missing thinning generally leaves too many; cannot raise an already deficient high-z bright parent |
| Incorrect/excessive upstream thinning | Removes too many | Possible in principle; no such BGS n(z) thinning established in inspected path |
| Fit ntilde(z), generate randoms or set weights | Changes expectation/reference/weighting | Does not regenerate the actual tracer population |

Previous common-sky high-shell raw r<19.5 parent counts are30,124SGC/65,661NGC,
versus41,243/94,985 successful real galaxies. Consequently, a selection-only
thinning step cannot turn that bright parent into the observed population.
The full raw r<20.2 pool is larger: changing the *physically justified* targeting
mapping could promote objects across the bright cut without rerunning the
underlying N-body simulation. That is a different operation from adding n(z).
Absence of a thinning stage alone is not a demonstrated cause of the shortage.

## Recommendation and next executable experiment

Do not regenerate the full production suite with an arbitrary supplied DESI
n(z) yet. Prepare one documented, already exposed phase as a reproducibility
experiment after resolving the column definitions. Reuse existing N-body/halo
products; a new gravity simulation is not indicated.

1. Obtain/identify the actual generator entrypoint, LF/HOD/colour tables,
   magnitude evolution definitions, random seed and any n(z) selection file.
   State whether n(z) describes intrinsic galaxies, targets or successful spectra,
   along with footprint and redshift measure. Pin the -.8 transformation.
2. Reproduce the existing phase first. Compare joined per-galaxy M/m/colour and
   shell counts before any changes, keeping source/config hashes and run logs.
3. If magnitude/targeting mapping is wrong, rebuild that layer from the available
   parent. If the physical population is insufficient, regenerate the galaxy
   population with a justified LF/HOD model. Apply each observation stage once.
4. Validate counts, magnitude/colour/size distributions, clustering and
   conditional environment/posterior coverage. Agreement of marginal n(z)
   alone is not an independent validation after it was used as an input.
5. If training observations change, regenerate fields and OOF summaries and
   refit/revalidate the posterior; assess encoder retraining with matched tests.
   Preserve independent phases and the current provisional VAC.

This is a recommendation, not a launched regeneration or repair. The current
request's diagnostic-only boundary remains in force.

## Reproducibility

workflows/p12a_vac/nz_magnitude_hypothesis.py uses only the saved512-row paired
SED sample and source snapshots, so it needs no new Slurm allocation. JSON
records source/sample/script hashes and the explicit conditional assumption.
The prior VAC implementation and audit artifacts are being committed along
with this report; unrelated E2E work is excluded.

## Authoritative LSS Abacus chain: user clarification

User identifies desihub/LSS as the authoritative processing chain from AbacusSummit
catalogues. Inspected upstream commit d942b990860e016e7558c740966fea3e5ce7e7f3;
immutable sources and SHA256 are in LSS_UPSTREAM.json. This pins inspected code,
not the historical revision executed for our products.

- prepare_mocks_Y3_bright.py applies the apparent-r cut and bypasses the dark
  STATUS nz selection for BGS. Its downsampling default is n.
- prepare_mocks_Y1_bright.py instead uses IN_Y1 for BGS and defaults downsampling
  to y. Crucially, prepare_script_bright_Y3.sh actually invokes the Y1 script,
  explicitly overrides downsampling to n, and sets rbandcut 19.5. Filenames alone
  therefore do not establish the executed footprint or defaults.
- ab2ndgen_bgsbright02loa_interactive.sh processes AbacusSummitBGS_v2 with
  mkCat_amtl.py and specdata loa-v1. It first creates plain BGS_BRIGHT full
  products, then constructs BGS_BRIGHT-02 clustering products.
- In mkCat_amtl.py, --nz y calls common.mknz and common.addnbar: measure weighted
  counts divided by shell volume, write an n(z) table, and attach density/FKP
  information. This is not an observed-n(z) input or galaxy-count matching step.
- The BGS_BRIGHT-02 branch with absmagmd=redshiftdep reads two polynomial
  coefficient files BGS_ANY_zmagcut_a/b.dat, splits at z=0.3, adds 0.078 and
  selects R_MAG_ABS below that redshift-dependent threshold. The plain
  BGS_BRIGHT product must not be conflated with this additional sample.
  The coefficient origin and applicability to our historical catalogue are not
  established by this source inspection.

Next closure: pin historical Abacus preparation and mkCat command/commit, join
raw -> forFA -> full TARGETIDs, and distinguish plain BRIGHT from -02/-21.5
products and Y1/DA2 footprint choices. Existing raw-parent deficit cannot be
explained by omitting an output-only --nz y operation. Do not substitute Holi
or introduce a fitted n(z) correction. No pipeline was executed or repaired.

Validation: all 15 source snapshot SHA256 checks and diagnostic syntax passed.
The 13 lightweight P12-A unittest checks passed. Graphify update and global
refresh succeeded inside cosmic_env. Verbatim source snapshots retain original
whitespace; the pre-existing production script trailing blank line is retained
to preserve its recorded source hash. No scientific qualification is implied
by these software checks.
