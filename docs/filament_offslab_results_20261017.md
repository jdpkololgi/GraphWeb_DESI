# Frozen P12-A: off-slab galaxy geometry

## Verdict

**Partial projection effect; hidden sheets are not a sufficient general
explanation.** Adding off-slab neighbours measurably weakens linearity relative
to planarity, even after neighbour-count matching. Most measurable suspicious
targets nevertheless remain line-dominant under the tested 3D tracer statistic.
This is not evidence that their tidal class must be filament, nor proof of model
failure or domain shift. No retraining or catalogue change follows.

## Bound inputs and estimator

Verified parent180977 galaxies includes all4459 original A06 targets. Parent
SHA2567a60a47f7bd0107611ece6c1e0f7338769151784d18bd8abf080994c4c3905fb;
full VAC/model/source checks from preceding extraction remain bound via READY.
This run freshly verifies parent hash and exact original TARGETID/P_WEB matching,
and independently reproduces the old projected anisotropy to1e-10 on valid rows.
No original posterior or checkpoint modified. Saved arrays retain exact target IDs,
all cohort indicators, per-view shapes/counts and3D shapes/normals/neighbour counts.

Main cohort:297 originally projected-line-like targets with Pfil<.2 and wall
argmax (293 also have Pwall>=.5). Controls:185 corresponding void argmax and173
original projected-line-like filament argmax; all4459 targets measured too.
Class probabilities define cohorts, not geometric truth labels.

Use unweighted covariance of neighbouring observed redshift-space galaxy positions,
excluding the target by row index. For descending covariance eigenvalues a>=b>=c:
L=(a-b)/a, P=(b-c)/a, S=c/a. L+P+S=1; these are *shape descriptors*, not posterior
class probabilities. Sheet-dominant means P>L and P>S. A finite-sample, noisy
shape proxy cannot define a true tidal wall. Positions use Planck18 Mpc.

Compare identical3D balls with/without neighbours outside |height|<10Mpc.
Radii5.173,10.346,15.519Mpc; min counts6/8/12. Main radius=7/.6766Mpc, but a
hard-radius covariance is NOT the Gaussian-smoothed matter tidal operator and
is not the model's full receptive field. Model context was never limited to the
20Mpc display slab. No intervention on model inputs was made here.

Parent padding, catalogue radial limits and inherited BOUNDARY_MPC>radius are
checked for targets. These guard truncation, not unobserved galaxies or a freshly
measured exact angular/veto selection. Parent includes unsupported VAC tracers.
Input preparation source uses broader .1–.6 context versus .15–.55 VAC; local
balls are checked inside VAC radial bounds. No claim of full upstream-context
ID replay, weighting correction, real-space deprojection or mass matching.

## Main results (minimum8 neighbours)

| Radius Mpc | Measurable of297 | Line-dominant | Sheet-dominant | Isotropic-dominant |
|---|---:|---:|---:|---:|
|5.173|20|16|3|1|
|10.346|244|191|49|4|
|15.519|297|211|59|27|

At main radius,53 targets are too sparse, not classified as non-sheets.
49/244=20.1% sheet-dominant; conditional95% spatial-block interval14.5–26.0%.
At15.519Mpc59/297=19.9%, interval13.7–26.3%. Changing main-radius neighbour
threshold from6 to12 gives55/278=19.8% and35/167=21.0%. Half-radius results
are sparse (zero targets at minimum12) and do not establish a small-scale result.

At main radius193 wall-cohort targets have>=8 neighbours in BOTH full and slab
balls. Mean paired change in P-L is+.1334, conditional95%CI[.0764,.1939].
Randomly thinning each full ball to its slab neighbour count (32 replicates,
seed17017) gives+.1154,[.0665,.1646]. This supports a genuine truncation effect
beyond raw neighbour-count increase, but does not remove finite-sample bias.
Paired transitions:27 of178 slab-line-dominant targets become sheet-dominant;
149 stay line-dominant and2 become isotropic. Of14 slab-sheet targets,12 remain
sheet and2 become line. The last target stays isotropic. Thus paired full-sheet
count39/193 versus14/193 slab-sheet: improvement, not majority conversion.

Among49 full-sheet cases,median |normal dot slab-normal|=.371 and32/49<.5,
consistent with some sheets intersecting the slab more nearly edge-on. This is
conditional on noisy covariance sheets; not an independently calibrated normal.

Full-ball filament-argmax controls are47/165=28.5% sheet-dominant at main radius,
and70/173=40.5% at1.5R. Their median main-radius count22 exceeds wall-cohort13
(void9). These unmatched controls show that sheet-like galaxy covariance is not
specific to model-wall labels; they are not a conditional model-performance test.

Intervals use512 paired spatial-block bootstrap draws on30Mpc Cartesian blocks,
conditional on fixed selected targets and proxy definition. They exclude model,
selection, mass/colour, redshift-space and between-realisation systematics;
nearest-neighbour volumes overlap and block adequacy is not independently proven.

## Thick slabs and orientations

Parent independently reselected for full depths5/10/20/40/80Mpc, y rotations
0/30/60/90deg and normal offsets-5/0/+5Mpc. Every original target remains tracked;
view membership/count losses are saved for all60 views. No favourable crop search.

Only3 wall-cohort targets are measurable across all five unrotated centred depths;
that registered all-depth comparison is not representative. Secondary paired
20/40/80Mpc comparison retains all297 wall targets,185void and173filament controls:

| Cohort | Line-like at20 | at40 | at80 |
|---|---:|---:|---:|
|low-Pfil wall argmax, n297|288|63|9|
|low-Pfil void argmax, n185|185|46|4|
|filament argmax, n173|161|41|8|

20Mpc wall count288 rather than original297 reflects lifting original RA/redshift
truncation when selecting the padded parent. The original statistic on exactly
original rows reproduces separately. Falling projected elongation in thicker
slabs is substantial but occurs in all three cohorts; superposition is not proof
of sheets. Rotation also changes target inclusion: at20Mpc/zero offset the wall
cohort valid counts297/30/17/10, line counts288/8/6/5 for0/30/60/90deg. Do not
interpret these unpaired rotation fractions as morphology changes of297 galaxies.

## Figures and provenance

- [Geometry and count control](figures/filament_offslab_20261017/01_geometry_comparison.png)
- [Thickness with fixed targets](figures/filament_offslab_20261017/02_thickness.png)
- [Same neighbours in perpendicular views](figures/filament_offslab_20261017/03_neighbour_views.png)

PDF versions adjacent. Three illustrated wall-cohort targets chosen by sorted
TARGETID quartiles25/50/75 among main-radius eligible targets, independent of
shape; IDs/probabilities and exact selection in EXAMPLES.json. Orange neighbours
are inside original slab; green crosses are off-slab. All panels equal physical
scales, no smoothing/jitter. Figures visually inspected for clipping/labels.

Scripts: workflows/p12a_vac/filament_geometry/{compare_offslab,validate_offslab,
plot_offslab,test_offslab}.py. Developer locks and launch commands:
`p12a_vac_operations.md`. Three synthetic shape tests pass: line/plane/sphere,
rotation+translation invariance, and edge-on slice of a true plane.183 independent
saved-result checks pass. Empty-bin validator bug fixed; no numerical rerun or
scientific values changed. Analysis59904230 completed0:0 in34s step; validator
59904267 failed on empty bin, corrected59904281 completed0:0. All allocations
released. Last small figure re-render uses saved results only.

Evidence JSON and1.9MB geometry arrays preserved under
`docs/evidence/filament_offslab_20261017/`; figure hashes adjacent. Full parent
stays in isolated Scratch. No fresh mock-phase access, inference or retraining.

## Decision

Do not explain away every low-Pfil string as an off-slab wall. There is measurable
projection sensitivity and a minority of sheet-like neighbourhoods, but the broad
visual-versus-tidal-class discrepancy is unresolved. Equally, do not call remaining
line-like galaxy neighbourhoods model errors. Tidal morphology depends on the
smoothed matter potential, not just local bright-tracer elongation; RSD, selection,
shot noise and target definitions remain distinct alternatives.

The next discriminating test is the SAME geometry statistic on matched mocks with
bound matter-derived tidal truth, separating how often true walls look like strings
from excess errors of frozen P12-A. TARP success under its tested mock distribution
does not decide either visual morphology or domain transfer. No retraining warranted
from this geometric screen alone; provisional VAC decision unchanged.
