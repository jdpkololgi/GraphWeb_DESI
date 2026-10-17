# True tidal walls can look like strings in the ph006 mock

The frozen mock test supports a mismatch between tracer appearance and tidal
class as the main explanation **in this development mock**. Most string-like,
low-Pfil wall predictions are true tidal walls, not misclassified filaments.
It does not establish the error rate of the corresponding DESI galaxies.

## Main result

Use R = 7/0.6766 = 10.346 Mpc, at least eight full-ball neighbours, and the
frozen projected anisotropy cut 0.5847730874520835. Fractions below restore the
saved shell sampling weights; brackets are conditional 95% cap/superblock
bootstrap intervals. Counts are unweighted.

| Cohort | Measurable rows | 3D line-dominated | 3D sheet-dominated |
|---|---:|---:|---:|
| True tidal walls | 4,521 | 67.8% [66.7,69.1] | 25.0% [23.8,26.3] |
| Projected-string true walls | 968 | 80.3% [78.1,82.8] | 19.0% [16.6,21.2] |
| Projected-string, Pfil<0.2, wall-argmax predictions | 772 | 80.6% [78.0,83.4] | 19.0% [16.3,21.7] |

Of those 772 predicted-wall strings, **643 are true walls, 28 true filaments,
101 true voids and none true clusters**. The weighted composition is 83.3%
wall [80.5,86.1], 3.6% filament [2.7,4.6], and 13.1% void. Their mean predicted
P(wall) is 81.9%. Actual wall frequency minus mean probability is +1.4 percentage
points, interval [-0.6,+3.6]; this does not show wall overconfidence in this
selected cohort. This is a conditional mean check, not a full calibration proof.

For comparison, all measurable low-Pfil wall predictions have 83.1% true walls
and 5.3% true filaments. The string subset does not show a descriptive excess of
filament contamination. These are overlapping, unmatched groups; no independent
or causal significance claim is made from that comparison.

Among 5,160 true walls with valid projected shapes and boundary clearance,
23.3% [21.7,24.8] pass the frozen projected-string cut. Thus true tidal walls can
look string-like in projection, and many remain line-dominated when all nearby
3D tracers are included. Covariance dominance is a descriptive tracer statistic,
not a tidal class or a visual classification by a human.

## Sensitivity and selection limits

At the main radius, minimum 6/8/12 neighbours gives 69.2/67.8/65.3% line-dominated
true walls and 24.3/25.0/25.7% sheets. For predicted-wall strings the true-wall
fractions are 82.5/83.3/86.8%, with filament fractions 3.7/3.6/4.2%.

At 1.5R, all 1,012 boundary-cleared predicted-wall strings have at least eight
neighbours: 82.8% are true walls and 3.7% true filaments. Their full geometry is
65.2% line-dominated and 25.9% sheet-dominated. The result therefore persists when
the main-radius sparsity exclusion is removed. Half-radius samples are much
smaller (35 predicted-wall strings at minimum eight) and are not decisive.

The original A06 sky/position cuts contain only 52 saved evaluation targets;
just four predicted-wall strings meet the main 3D criteria. All four are true
walls, but there are only two spatial blocks. Those counts are illustrative;
the numerical bootstrap intervals stored for this tiny cohort are not reliable
uncertainty estimates and should not support conclusions.

The headline results instead use **target-centred, globally oriented 20 Mpc
slabs** across both mock caps at 0.20 <= z < 0.30. This retains the radius,
statistic and thresholds, but is not an exact replication of A06 selection.
Target centring makes off-slab losses small at the main radius; these results
test the relationship of geometry to truth, not the magnitude of the original
projection effect. Original-slab neighbours use a padded parent, rather than
clipping at the displayed rectangle. True-wall versus true-filament geometry
fractions are not density/count matched and should not be interpreted as intrinsic
morphological differences between those populations.

## Provenance and execution

Only already-exposed ph006 evaluation products were accessed. Its 50,000 saved
rows contain 21,052 targets at 0.20 <= z < 0.30; 9,010 carry true-wall labels.
Neighbour searches use 4,968,208 successful-redshift canonical context tracers
at 0.15 <= z < 0.55, including tracers with no valid training truth. No new model
inference, retraining, or sealed phase was used. This phase previously supported
model selection, so this is a development diagnostic, not blind validation.

Checks bind the dataset and checkpoint hashes, 512 frozen posterior draws per
row, evaluation folds, canonical parent-ID identity, position/index hashes, and
canonical truth-array hash. Saved sample eigenvalues match the canonical truth
at those parent IDs. Ordered-increment inversion matches the saved eigenvalues;
class counts use strict lambda>0.2 per posterior draw. Observed Cartesian
coordinates match Planck18 distances to the specified 0.15 Mpc tolerance.

The existing target contract specifies native CACTUS matter tidal labels with
7 Mpc/h smoothing, ordered eigenvalues and threshold 0.2. These are attached
native host-environment labels, not a tidal operator evaluated on observed
redshift-space tracers. This audit verifies stored-target correspondence; it
does not recompute the underlying matter tensor or revalidate every native
annotation operation. Redshift distortions, sparse tracer sampling and the
non-equivalence of covariance and tidal operators remain relevant.

Shape radii are 0.5/1/1.5R with neighbour minima 6/8/12. Self is excluded by
parent ID. Conservative random-support boundary clearance exceeds the projected
query bounding-sphere radius; it is not an exact angular-mask completeness proof.
Full/slab and 32 count-matched full-ball subsets are retained in the evidence.
Intervals use 512 cap/superblock draws, conditional on one phase and frozen model;
they do not cover simulation, selection or DESI domain-shift systematics.

User-approved CPU allocation **59904645**, nid200009, desi/interactive,
64 logical CPUs, 30-minute limit, Scratch license. Completed and released after
3m27s; analysis step 1m37s and validator/figure step 3s, both exit 0:0.
Three synthetic geometry tests and **1,080 independent saved-array checks pass**.
Figure visually checked; unrelated jobs untouched.

- Runner: `workflows/p12a_vac/filament_geometry/mock_truth.py`
- Validator/figure: `workflows/p12a_vac/filament_geometry/validate_mock.py`
- Evidence: `docs/evidence/filament_mock_truth_20261017/`
- Figure: `docs/figures/filament_mock_truth_20261017/mock_truth.png` (also PDF)
- Scratch: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_mock_truth_20261017_v1/`

No retraining is justified by the string appearance alone. The evidence favours
true-wall tracer strings in this mock; whether that explanation transfers to
Loa remains a separate question.
