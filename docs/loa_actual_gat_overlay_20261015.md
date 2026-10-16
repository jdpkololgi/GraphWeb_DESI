# Actual legacy GAT versus Abacus on the same Loa galaxies

The four newly fitted **actual GAT classifiers** do not recover a dramatically
stronger property association than Abacus on the tested Loa samples. This
holds using the historical observed-feature normalization and the recovered
mass-completeness method. The cuts mainly change which stellar-mass population
is being described. Frozen-transform failures make large raw contrasts
unreliable evidence of successful transfer.

This result is from one bounded fit per arm, not a convergence study, a
multi-seed architecture ranking, or a reanalysis at the old z<=0.06 selection.
The current windows are 0.20–0.30, 0.20–0.25 and 0.25–0.30. The earlier MLP
selection experiment is separately documented in
`loa_selection_sensitivity_20261015.md`; its trace-quartile contrasts are a
different estimand and must not be substituted for the class contrasts here.

## Models, identity and verification

Used the canonical Illustris `SimpleGAT(10,4,num_heads=4)`: four GATv2 layers,
the original residual/LayerNorm/dropout structure, edge length and ten graph
features. Rebuilt normal and thinned alpha-complex graphs and all graph
features. R2/R7 labels come from the same verified TNG300-1 matter truth,
threshold 0.2. The sparse population retains 37,545 of 283,969 tracers.

The normal graph exactly reproduces all **3,890,098 directed edges** of the
archived GAT training cache. Positions agree to 1.71e-13 Mpc; all ten
full-fit transformed features agree to 5.66e-7. Independent NetworkX weighted
clustering and direct-centered covariance checks pass for each graph.
Trained destination-block inference agrees with the canonical full forward;
Loa layerwise versus exact four-hop inference differs by at most 3.82e-6 in
logits. Input hashes, held-out IDs, labels and saved predictions pass replay.
This establishes actual model/input identity, not observational calibration.

Training used one paired seed, 1,200 Adam minibatch updates with 16 destination
nodes and their complete four-hop context. Train-only scaling/class weights,
spatial label splits, and validation-loss checkpoint selection. No sampled
neighbor truncation. The 2,048 held-out galaxies are identical across arms:

| GAT | Best update | Held-out balanced accuracy | Macro F1 |
|---|---:|---:|---:|
| R2 normal | 1000 | 71.03% | 0.648 |
| R2 thinned | 1100 | 61.25% | 0.556 |
| R7 normal | 900 | 67.56% | 0.628 |
| R7 thinned | 900 | 71.36% | 0.660 |

The R7 ordering does not establish a benefit of thinning: optimization/seed
uncertainty is unresolved. Graph contexts are transductive; sparse training
has fewer labelled targets. The historical nonperiodic graph is retained.
Full contract: Illustris `docs/legacy_gat_density_20261015.md`.

## Observational comparison

All models infer on the same 103,788 Loa query IDs using one graph of 587,896
padded observed tracers, with 3,052,109 undirected edges. No Loa thinning.
Historical DESI alpha=6.90626 Mpc is retained for all four models; sparse TNG
training uses its density-dependent alpha=13.55655 Mpc. This inherited graph-
scale mismatch is a remaining transfer factor, not an isolated density test.
Outcome cuts do not change the graph, predictions or normalization.

The valid CG15/magnitude-limited parent has 90,952 galaxies. The historical
faintest-20% / 90th-percentile limiting-mass method retains 52,620. Estimate
the limit before applying the outcome mass range. Within dz=0.01, use
logM_lim=logM+0.4(r-r_limit), with extinction-corrected BRIGHT limits 19.5
South/19.54 North. This is not a correction for all survey selection effects.

The primary historical convention fits Box–Cox to this observed graph's
features and recenters them. It reproduces the old default **operation**,
not the old catalogue's moments. The sensitivity convention applies the
unchanged training transform. Neither uses galaxy-property outcomes.

Soft curves use GAT softmax scores and Abacus P_WEB weights; GAT class scores
are not calibrated posterior probabilities. Hard curves use argmax for
every model, including Abacus. Mean-eigenvalue classes from the earlier MLP
comparison are not mixed into this estimator.

## Completeness does not restore a large missing contrast

Probability-weighted cluster-minus-void quenched contrast, percentage points,
at a **single shared parent/cut stellar-mass and redshift mixture**:

| Model | Parent | 90% completeness cut |
|---|---:|---:|
| Abacus P12-A | 9.460 | 9.647 |
| GAT R2 normal | 8.732 | 8.590 |
| GAT R2 thinned | 5.893 | 5.894 |
| GAT R7 normal | 6.031 | 6.015 |
| GAT R7 thinned | 7.595 | 7.644 |

This reference uses 67 common 0.1-dex mass / 0.01-z cells containing 52,340
parent rows; selected weights are zero for cut failures. Minimum Kish count
20 per group/cell. Uncertainty uses 128 paired draws of 21 sky blocks, with
16–84% intervals in RESULTS.json. These are conditional sampling intervals,
not model/SED/transfer uncertainty. Small differences do not establish an
architecture ranking.

If each selected sample gets its own common reference mixture, the apparent
before/after values instead become:

| Model | Uncut mixture | Cut mixture |
|---|---:|---:|
| Abacus | 7.621 | 9.145 |
| GAT R2 normal | 6.810 | 8.447 |
| GAT R2 thinned | 4.867 | 5.865 |
| GAT R7 normal | 4.688 | 5.772 |
| GAT R7 thinned | 6.257 | 7.420 |

These panels contain 83,439 and 48,772 common-cell rows. The increase is not
an estimate of the completeness cut's effect at fixed population. The lower
redshift half shows somewhat larger R2 contrasts, but does not restore a
large old-GAT advantage: its cut-mixture Abacus/R2-normal/R2-thin values are
9.89/9.33/6.44 pp; the upper-half values are 8.27/7.42/5.20 pp.

Raw soft contrasts actually decrease with the mass cut: Abacus 16.59→12.56,
R2 normal 12.48→6.81, R2 thin 11.04→6.06, R7 normal 6.54→3.54,
R7 thin 9.06→4.70 pp. All seven property curves are preserved in the PDF.

## Frozen-transform failures and boundaries

With the frozen training transform, normal R2 places 82.05% of selected
galaxies in voids and none in clusters. Normal R7 places 95.45% in voids
and reverses the raw quenching association. Sparse R2 places 65.04% in
clusters but only 0.046% in voids; its superficially strong 24.24-pp raw hard
contrast relies on that very small void group. Sparse R7 places 92.93% in
walls. There is no shared four-class controlled support across the five
models, even for soft scores, in any tested window/cut. Those plots remain
explicitly unavailable; no controlled frozen-versus-historical contrast is
claimed. Historical normalization changes class occupancy, not truth labels.

The common parent retained the preceding local-feature geometry mask, which
is not by itself a full GAT-context guarantee. A conservative additional
margin >12alpha=82.875 Mpc covers four message hops plus two feature-context
hops. It retains 40,571 galaxies, of which 31,260 share controlled support.
Historical-normalization contrasts are Abacus 8.68, R2 normal/thin 7.79/5.26,
R7 normal/thin 5.15/6.94 pp. The main conclusion persists on this interior
sample. This does not remove RSD, targeting or large-scale survey selection.

## Artifacts and execution

Source: `workflows/presentation/loa_actual_gat_overlay_20261015.py`.
Review bundle: `docs/figures/loa_actual_gat_overlay_20261015/`.
Full outputs: `~/fmsu_tokyo_20261010/loa_actual_gat_overlay_20261015/`.
The 32-page PDF contains seven-property soft raw/controlled overlays plus
hard-label quenching checks, for both transforms and all three redshift
windows. CSV, complete JSON curves/support/occupancy, validation and hashes
make empty groups and negative results inspectable.

All stages completed sequentially using the explicitly authorized one-thread
login fallback. No Slurm job or production VAC was changed. `graphify` is
unavailable, so its index could not be refreshed; direct source/receipt
verification was used. No model was tuned using observed trend strength.
