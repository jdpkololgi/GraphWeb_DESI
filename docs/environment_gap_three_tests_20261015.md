# Why did the galaxy–environment relation change?

**Evidence synthesis • 15 October 2026 • TNG GAT versus Abacus P12-A**

## Answer and interpretation of the Gemini discussion

I agree with the broad warning: a stronger property–inferred-environment association does not establish that the old GAT learned the physics better. We changed the galaxy population, physical environment target, observation process, feature normalization and estimator. The evidence now supports **population composition and target scale as substantial mechanisms**, and identifies **transfer/normalization failures**. It does not identify an additive causal budget for the historical observed gap.

Several stronger statements in the quoted discussion are not justified. Not every change necessarily flattens every statistic. A 7 Mpc/h filter does not remove all halo-related information; mass and environment remain correlated. A GAT failure at high redshift does not uniquely identify density dependence: normalization, graph scale, selection, RSD and training coverage all change. Geometry predicting density is not itself illicit property leakage; the tested models do not ingest stellar mass, colour or SFR. Mock-conditional calibration is not established real-survey calibration, and posterior means or probability weighting do not obey a universal attenuation formula.

The new tests are particularly informative: **the high-mass TNG population has much weaker true property contrasts, but selecting the narrowest observed posteriors weakens—not restores—the signal.** These findings support the selection hypothesis while contradicting the claimed automatic benefit of certainty cuts.

The discussion's satellite fractions were estimates, not measurements of our samples. In the held-out TNG population, the measured satellite fraction is **40.62%**, falling to **33.75%** above log M=10.5. These are z=0 simulation fractions; the DESI satellite fraction remains unmeasured here. A high-mass cut does not turn this into a sample dominated almost entirely by centrals. Nor was the old observational selection exactly physical M>10^9: the recovered FastSpecFit convention implies approximately log M≥9.305 after its h conversion.

Published work supports distinct mass and environmental associations and an important satellite contribution, but does not specify the amplitudes or satellite fractions of these particular samples: [Peng et al. 2010](https://arxiv.org/abs/1003.4747), [Wetzel et al. 2013](https://arxiv.org/abs/1206.3571), [Kovač et al. 2014](https://academic.oup.com/mnras/article/438/1/717/1040681). The conclusions below are based on our measured controls, not on importing those studies' numerical effects.

## Reading the measurements

Void, wall, filament and cluster mean respectively zero, one, two or three tidal eigenvalues above threshold 0.2. Eigenvalues are ordered λ₁≤λ₂≤λ₃; trace is their sum. R2/R7 denote Gaussian matter-field smoothing radii of 2/7 h⁻¹ Mpc. A count aperture labelled N7 is a different operation from an R7 tidal field.

Quenching means sSFR<10⁻¹¹ yr⁻¹. A contrast of 10 percentage points (pp) means a 0.10 difference in quenched fraction. **Class contrasts** compare cluster with void; **quartile contrasts** compare high with low trace within defined strata. Soft class curves use membership weights; hard curves use labels. These estimands are not interchangeable. GAT softmax scores are not certified calibrated posteriors.

Every numerical comparison below identifies its population and control mixture. A change of common support changes the question. Bootstrap intervals describe spatial sampling conditional on the existing fits; they do not include training uncertainty, SED systematics, simulation uncertainty or survey-domain calibration.

## Test 1 — True TNG environment: scale and mass selection before ML

**Status: complete for the proposed high-mass-cut alternative.** No mock flux-limited lightcone was built. The new calculation reuses the verified common R2/R7 matter field and 66,691 held-out TNG300-1 z=0 galaxies with 9≤log M<11.5. There are 9,338 above log M=10.5. Query properties, positions and truth remain fixed. Quenched fraction, g−r and capped log sSFR are measured directly, with no predicted environment.

Trace quartiles are assigned within fixed 0.25-dex mass bins. Each contrast standardizes its high/low groups to that population's mass-bin proportions; the high-mass sample uses the corresponding subset of the same bins and quartiles. This tests changing the population, not an effect on an individual galaxy. The cutoff is the user's suggested proxy; it is **not** an equivalent r<19.5 selection, and has neither z≈0.3 evolution nor DESI selection effects.

| Population | N | R2 contrast [95% interval] (pp) | R7 contrast [95% interval] (pp) |
|---|---:|---:|---:|
| Full: 9≤log M<11.5 | 66,691 | 49.50 [48.03, 50.96] | 32.73 [30.48, 35.92] |
| High mass: 10.5≤log M<11.5 | 9,338 | 12.03 [9.90, 14.11] | 10.11 [7.68, 12.59] |

Smoothing alone reduces the full-population quartile contrast by 33.88%. The high-mass cut reduces it by 75.69% at R2 and 69.12% at R7. Within the high-mass population, changing R2 to R7 reduces the contrast by about 16.02%. These conditional effects interact; do not add them, and do not claim that they explain those percentages of the observed old/current gap.

The earlier same-support physical-class comparison was 52.04/47.45 pp at R2/R7, a much smaller scale difference than the quartile comparison. The new supplementary physical-class results in test1.json use each arm's own eligible cells and therefore must not be subtracted as a fixed-support scale effect. This is direct evidence against treating a single reduction as universal.

**Satellite census:** a central is a selected subhalo whose ID appears in GroupFirstSub in the complete TNG snapshot-99 group catalogue; others are satellites. The measured 40.62→33.75% change is consistent with some population shift, but does not isolate satellites as the cause of the contrast reduction. A within-mass central/satellite decomposition and matched observational groups would be needed for that attribution.

**Existing results assigned here:** native R2 operator recovery (effective R=2.02662), 8,192/8,192 sampled field-label matches, two discrepancies among 283,969 cached labels, mesh/interpolation checks, the same-field R2/R7 property experiment, and the new high-mass/census diagnostic. Full receipts and the original methods report are included in the companion document's expandable evidence sections.

<!-- FIG_TEST1 -->

## Test 2 — Both model families on the same observed galaxies

**Status: the higher-redshift alternative is completed on 0.20–0.30 and its two halves.** It is not a full 0.15–0.50 GAT run, nor a 0.15–0.20 paired run. The existing test is already sufficient to ask whether a TNG-trained GAT necessarily produces the large historical contrast on a current Loa population. It is not a matched-truth comparison of production architectures.

Four actual legacy-architecture GATs were trained at R2/R7 and normal/thinned TNG density. The graph/features and canonical four-layer SimpleGAT identity were verified, including exact replay of 3,890,098 directed training edges. Each fit used one seed and 1,200 updates; this is a bounded comparison, not demonstrated convergence. Held-out balanced accuracies were 71.03%, 61.25%, 67.56%, 71.36% for R2 normal/thin and R7 normal/thin respectively. They cannot establish a beneficial effect of thinning at R7.

All models were evaluated on identical Loa query IDs using the same observed tracer graph. On a single common parent/completeness-cut mass–redshift mixture, the soft cluster-minus-void quenched contrasts are:

| Model | Parent (pp) | 90% limiting-mass cut (pp) |
|---|---:|---:|
| Abacus P12-A | 9.460 | 9.647 |
| GAT R2 normal | 8.732 | 8.590 |
| GAT R2 thinned | 5.893 | 5.894 |
| GAT R7 normal | 6.031 | 6.015 |
| GAT R7 thinned | 7.595 | 7.644 |

There is **no recovery of a dramatic old-GAT advantage** in this test. A conservative >12α graph-context interior cut preserves that conclusion (Abacus 8.68; GAT R2 normal/thin 7.79/5.26; R7 normal/thin 5.15/6.94 pp). All seven property panels, both redshift halves, hard/soft choices and both normalization conventions are preserved in the 32-page embedded supplement.

The historical normalization operation fits Box–Cox to the observed graph and recenters it. Freezing the training transform instead causes extreme occupancy failures: normal R2 assigns 82.05% to voids and no clusters, normal R7 95.45% to voids, thin R2 65.04% to clusters but only 0.046% to voids, and thin R7 92.93% to walls. No five-model four-class common controlled support survives. A large raw contrast from a tiny void group is not a successful scientific recovery.

This failure identifies transfer sensitivity, not density alone. The historical Loa graph uses α=6.90626 Mpc, whereas the thinned TNG training graph uses α=13.55655 Mpc. Label scales, graph scales, normalization and survey observation processes still differ. The newly trained GATs are not the exact recovered historical weight file; that file was separately replayed to identify the saved old catalogue.

The earlier geometry-only MLP transfer experiment is retained separately: it is not the actual GAT. Its rank contrasts also show that changing the population reference can create an apparent completeness-cut improvement. At a fixed parent/cut mixture the Abacus/MLP contrasts barely change. Its full 21-page selection supplement and initial overlay figures are included.

**Existing results assigned here:** checkpoint identification/replay; actual GAT architecture/input parity; four GAT training arms; same-ID Loa inference; normalization stress test; boundary sensitivity; three redshift windows; seven-property curves; initial MLP transfer and physical/rank/soft-Abacus estimator comparisons. The original saved VACs themselves have no usable same-galaxy overlap: seven shared IDs, zero consistent redshifts.

<!-- FIG_TEST2 -->

## Test 3 — Mass/redshift controls and actual posterior uncertainty

**Status: completed now across 0.15–0.50 in the established angular wedge, including the requested 0.15–0.25 focus.** This expands the redshift range of the earlier property diagnostics. It is not an all-sky analysis. The new extraction contains 236,774 SUPPORTED rows at RA[120,160), Dec[14.5,30.6). The quenching-only sample retains 221,002 unique, unambiguous matched GALAXY rows with finite positive CG15 mass, nonnegative SFR, 9≤log M<12, valid eigenvalues and normalized P_WEB. It does not inherit the seven-property AGE/AV/FAGN requirements, BRIGHT re-selection or the earlier local-aperture geometry mask; therefore its amplitudes need not exactly match earlier 0.20–0.30 tables.

Controls use **0.2 dex in stellar mass and 0.05 in redshift**. The available VAC stores eigenvalue quantiles rather than standard deviations. We therefore use **w₁=(Q84(λ₁)−Q16(λ₁))/2**, the central-68%-interval half-width, explicitly not σ unless the posterior is approximately Gaussian. Select exactly the lowest floor(N/2) widths inside each mass/redshift cell, resolving ties by TARGETID. No galaxy-property outcome determines the selection.

For the soft-class contrasts, parent and selected samples share one parent reference mixture over cells with Kish effective count≥20 in every class in both samples. Width selection thus cannot improve the comparison merely by changing the mass/redshift cell weights. There are 200 paired spatial block draws; listed intervals are 95% conditional intervals. Boundary sensitivity additionally requires the stored boundary distance>20 Mpc. This coarse stored distance is a sensitivity check, not a proof that every model context is unaffected by masks or collisions.

| z range | Common parent N | All eligible (pp) | Low-width half (pp) |
|---|---:|---:|---:|
| 0.15–0.50 | 218,556 | 7.45 [6.47, 8.20] | 1.62 [1.12, 2.07] |
| 0.15–0.25 | 115,518 | 7.68 [6.27, 8.98] | 1.05 [0.43, 1.57] |
| 0.15–0.20 | 61,262 | 7.76 [6.05, 9.76] | 0.79 [-0.12, 1.76] |
| 0.20–0.25 | 54,256 | 7.58 [5.86, 9.17] | 1.33 [0.37, 2.24] |
| 0.25–0.30 | 47,181 | 7.65 [5.83, 9.35] | 2.12 [0.77, 3.18] |
| 0.30–0.35 | 32,647 | 7.40 [5.85, 8.87] | 2.91 [1.71, 4.16] |
| 0.35–0.40 | 15,667 | 6.78 [4.53, 9.19] | 2.15 [0.22, 4.14] |
| 0.40–0.45 | 5,927 | 4.24 [1.42, 7.37] | 0.02 [-2.97, 3.08] |
| 0.45–0.50 | 1,616 | 4.61 [-1.25, 10.33] | 2.63 [-3.87, 8.71] |

The low-width cut **reduces** the property contrast. In the 0.15–0.25 focus it falls 7.68→1.05 pp; with the boundary cut it is 7.73→1.01 pp. The full-range result is 7.45→1.62 pp. Thus the proposed “remove uncertain nodes and recover the missing relation” explanation is not supported by this executed diagnostic.

Continuous checks agree in direction: rank-transform λ₁, trace, quenched status, capped log sSFR, mass and z; subtract cell means; then regress out the residual ranked mass/z terms before correlating. In 0.15–0.25, partial trace–log-sSFR correlation changes from −0.0757 to −0.0403, and λ₁–log-sSFR from −0.0688 to −0.0313. These are descriptive partial rank correlations, not causal coefficients; the minimum-40-row cells are selected separately per arm, so use the shared-mixture class comparison for the primary before/after statement. High-z bins have much less support and wider intervals.

Selecting narrower posteriors is itself environment selection, not a calibration repair. Width can depend on the inferred density, class and prior; conditioning on it can restrict dynamic range or induce selection bias. The measured change is large: at 0.15–0.25, mean cluster weight falls from 5.50% to 0.185%, while void weight rises from 25.07% to 48.63%. The median inferred trace moves from 0.294 to −0.031. Thus the selected soft cluster curve largely relies on small tail weights, even when its Kish count passes; it does not represent a comparable population of confidently identified clusters. The saved selection profiles report these values before common-support restriction. This result does not demonstrate that wider posteriors are better, that the model is calibrated on DESI, or that uncertainty never attenuates a relation.

<!-- FIG_TEST3 -->

### Earlier selection, weighting and property controls

| Diagnostic | Existing result | What it establishes |
|---|---|---|
| Common CG15 property source, own controlled mass mixture | Old/current soft quenching 21.87/8.97 pp | A difference remains after consistent property measurement; redshifts differ. |
| Shared between-sample mass reference | Soft 14.37/8.81 pp | Much of the initial gap is sensitive to stellar-mass composition/support. |
| Hard labels on that shared reference | 14.89/11.85 pp | Probability weighting contributes to the displayed gap; not calibration evidence. |
| Main BRIGHT restriction in old sample | 21.872→21.925 pp | Target-program mix is small for this comparison. |
| Same angular wedge | Only 172 rows meet strict shared controlled support | No robust isolated sky contribution. |
| Current limiting-mass cut, paired mixture | 9.459→9.425 pp; +0.2 dex cut 8.940→8.562 | This completeness cut does not restore the missing contrast at fixed control mixture. |
| Maximum class probability≥0.7 | Paired current soft 9.735→11.119 pp | A different confidence statistic selects stronger trends; thresholds .85/.95 fail common support. Not the new eigenvalue-width test. |
| Sharpen probabilities into hard labels | Old 21.872→22.667; current 8.970→11.764 pp | Estimator sensitivity on identical rows; no correction for uncertainty. |
| CG5 replacing CG15 | Seven trend directions persist; quenching 8.970→9.573 pp | This SED choice is not the dominant explanation on the tested sample. |
| Completeness rationale recovered | Faintest-20%/90th-percentile curve reaches log M=9 at z=.055–.06 | Supports the remembered motivation; exact historical decision receipt is missing. |

The often-quoted gap reduction from 12.90 to 3.04 pp combines mass-mixture/support and hard/soft-estimator changes. It is **not** a causal decomposition. The low-z and current catalogues have disjoint redshift support, so no weighting scheme can turn them into a same-redshift comparison.

## Additional test 4 — Density interventions, independent environments and provenance

These findings do not fit entirely inside the proposed three tests and are retained explicitly.

### 4A. Change the tracer information while holding truth fixed

TNG random thinning from 283,969 to 37,545 tracers, with fixed held-out galaxies/truth and retraining, gives mass-controlled quenched trace-quartile contrasts 48.52→45.02 pp at R2 and 33.75→32.64 pp at R7. The losses are 7.22% [5.49,8.97] and 3.28% [1.64,5.42], respectively, over three paired fits with spatial bootstrap. Frozen dense-trained regressors fail much more severely. This separates reduced information from failure to adapt and does not reproduce luminosity selection, RSD or fibre assignment.

Before training, the fixed-query TNG count experiment showed rho(N2,true R2 trace) .891→.723–.733 at 10% tracer retention; the zero-N2 fraction grows 7.6→45.8–46.3%. At 1% retention rho is .431–.447. All density levels (.3,.1,.03,.01), count apertures and three-seed ranges remain in the full ledger.

Abacus has its own fixed-query/truth count intervention. Across z=.15–.25/.25–.35/.35–.45/.45–.55, full-density N7/true-R7-trace correlations are .760/.708/.620/.414. At 10% retention they are .566–.579/.423–.434/.289–.302/.120–.137. Selection-adjusted N7 improves the full-density shell values to .791/.754/.677/.496. These demonstrate information loss, but are not numerical predictions of the observed property-gap amplitude.

The existing MT4 test adds actual Faint neighbours at fixed Bright targets/truth and compares angularly scrambled Faint context. Actual context beats scrambled context in every tested shell/rotation; macro λ₁ R² is .6074/.5394 versus .4753/.5187 for rotations 0/2. Bright-only is .5070/.5197. This is evidence for added spatial information, but the rotation-2 promotion safeguards failed. No production promotion follows.

### 4B. Independent observed environments and geometry

On 91,066 common observed rows, soft web classes give a 7.42-pp quenched contrast; independent N4/N7 and projected Σ5 give 7.82/7.68/9.52 pp. Inferred trace gives 6.35 pp. N1 gives 3.30 pp but 76.1% of all query galaxies have zero neighbours inside 1 h⁻¹ Mpc. This is a sparse indicator, not evidence that small-scale physics is intrinsically weaker.

Adding N7 control on the same 4,539 galaxies reduces the soft-class association from 4.56 to 2.17 pp. Compactness retains about 6.23 pp; the 3D quadrupole is −2.67 pp, while sky-plane quadrupole and inferred shear intervals include zero. This is not causal mediation: finite-width density control, geometry selection and redshift-space effects remain. The original point-density, multiscale, projected, concentration, quadrupole, shear and eigenvalue-distribution figures are all retained. Earlier 7.69/7.36-pp point-density/class results use a different support and are not contradictory revisions of the model.

### 4C. Identity, historical artifacts and low-z mock mapping

The saved old VAC passes the corrected north-then-south ordering audit; 487 repeated IDs were excluded from controlled comparisons. There are 193,706 unique consistent CG15 matches. Current inference/VAC/source crosschecks pass the recorded ID census and sampled draw/merge replays. The January 15, 2026 checkpoint reproduces the saved old probabilities to 1.19×10⁻⁷. **This identifies the available saved VAC, not necessarily the exact plot shown last year.** Historical plot/colour/mass reproductions are retained; the date/provenance distinction limits claims about the original slide.

Low-z Abacus BOX_INDEX=−1 is a separate population/indexing branch, not proof of a different epoch. Substituting snapshot .1 for .2 does not fix raw host matching (0.32% versus 0.37% in the probe); a high-z .2 control matches 92.7%. Reconstructing public halo filters matches only .53/.55% for raw/cleaned masses. Exact production crosswalk recovery remains incomplete. Do not use this unsupported low-z mapping to retrain an ostensibly matched R2 Abacus model and claim a resolved comparison.

### 4D. Alternatives and proposed work, kept distinct from completed results

The existing simulation audit covers Abacus base/high/highbase/small, Uchuu, Rosella, GLAM, MillenniumTNG, FLAMINGO and TNG300. The complete table, primary-source links and local support limitations are preserved in the original ledger. No verified candidate currently supplies every desired combination of independent phases, high-resolution galaxy physics, BGS survey selection and a complete galaxy-to-matter truth crosswalk. Alternative simulations and high-resolution pretraining remain options, not completed demonstrations of improvement.

## What can now be concluded, and what remains unidentified?

1. **Demonstrated:** changing the mass population and target smoothing can substantially change the true relation even without ML. The new high-mass TNG test is the strongest direct support for the physical/population explanation.
2. **Demonstrated:** the observed historical gap is sensitive to population reference and hard/soft estimation. Same-current-galaxy actual GAT comparisons do not restore a large old advantage. Normalization and graph-domain transfer can radically change class assignments.
3. **Demonstrated:** sparse tracers lose information, but adapted TNG predictors show relatively modest property-contrast losses under random thinning. Density alone has not reproduced the observed gap.
4. **Negative result:** neither the tested completeness cuts at fixed mixture nor the new low-posterior-width selection restores the large missing contrast. The latter weakens it substantially.
5. **Unidentified:** the separate contributions of real satellite loss, evolving galaxy physics, luminosity-dependent tracer bias, fibre assignment, RSD, survey masks and simulation/HOD differences. There is no observational matter truth here with which to certify the stronger relation as more accurate.

A defensible scientific statement is: **“The change in the measured relation is substantially affected by population and estimand choices; controlled TNG truth confirms strong mass-selection and scale dependence, and same-galaxy transfer tests expose normalization/domain sensitivity. The evidence does not support superior physical learning by the historical GAT, nor tracer density or posterior broadening as a sufficient single explanation.”**

To identify a unique dominant causal contribution to the exact historical slide, we would still need its exact artifact/selection and a validated paired galaxy/property/truth mock subjected sequentially to luminosity selection, RSD, masking and assignment, with both methods tested on identical targets. That is a new validated survey-forward experiment; the three diagnostics cannot supply that counterfactual merely by slicing existing predictions. No unsupported low-z truth, new production fit or reserved-phase access was used to fill this gap.

## Reproducibility and document coverage

The self-contained HTML companion includes this synthesis, new plots/tables, the complete earlier test ledger, original methods/result reports, all distinct PNG/SVG figures in the named investigation bundles, and embedded original multipage PDF supplements. Expandable source sections preserve historical wording and are clearly marked: older “open” statuses are superseded where the synthesis records a completed follow-up. The catalogue and model families remain separate.

New code: `GraphWeb_DESI/workflows/presentation/environment_gap_three_tests_20261015.py`. New result receipts: `GraphWeb_DESI/docs/figures/environment_gap_three_tests_20261015/`. The extracted observational cache stays outside the repository under `fmsu_tokyo_20261010/environment_gap_three_tests_20261015/`. Validation covers direct weighted-mean reconstructions, selected source-row parity, TNG agreement with the previous full-population result, image integrity and document inventory/hashes. Source FITS identity is bound by the recorded path/size/mtime, row census, fresh row checks and extracted hash; the full large FITS hash was not recomputed.

Execution used the user's authorized sequential one-thread login fallback. No Slurm job, production catalogue or scientific redshift policy changed. Graphify remains unavailable; direct scoped source/receipt inspection was used. The requested text-and-images alternative is supplied as a self-contained HTML file; the report-to-PDF skill's required exporter is unavailable, so no PDF conversion is claimed. Existing plot PDFs are embedded unchanged.
