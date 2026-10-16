# Shareable mass-selection and thinning figures; next causal tests

## Completed figures

Bundle: `/global/u2/d/dkololgi/fmsu_tokyo_20261010/supervisor_environment_figures_20261015/`.
Recommended: `tng_truth_prediction_quenched.pdf` and `loa_completeness_explanation.pdf`.
Also: TNG colour and log-sSFR curves; contrast summary including frozen-transfer
predictions; all seven raw Loa property curves. All have PNG counterparts.

TNG uses fixed held-out queries, matter truth, cached predictions, 0.25-dex
mass bins and fractional quartiles. Population cut: 9<=logM<11.5 versus
10.5<=logM<11.5. Mass cuts select queries only: they do not remove low-mass
neighbours or retrain the model. Random thinning instead changes context from
283,969 to37,545 tracers; predictions are from three paired geometry-only MLP
fits, not GAT/P12-A. 200 paired spatial bootstrap draws; fit uncertainty is
not included. 280 independent block/direct-mean checks pass; full-population
truth contrasts reproduce the original fractional-quartile analysis exactly.
Small differences from the latest threshold-based high-mass diagnostic
(12.03/10.11 versus12.04/10.10pp) reflect fractional rank boundaries.

| Population / scale | True contrast | Dense prediction | Sparse retrained | Frozen dense model on sparse inputs |
|---|---:|---:|---:|---:|
| Full / R2 |49.50|48.52|45.02|34.38|
| High mass / R2 |12.04|9.96|10.85|8.87|
| Full / R7 |32.73|33.75|32.64|29.13|
| High mass / R7 |10.10|9.70|9.68|8.56|

Entries are quenched-fraction Q4−Q1 differences in percentage points.
Do not equate these rank contrasts with observed cluster−void contrasts.

Loa figures replot verified existing P12-A class probabilities/CG15 summaries
on z=.20−.30, not new inference. The parent already requires mass>=1e9
(and<1e12), BRIGHT magnitude and seven-property/geometry eligibility.
The 90% limit rises from logM10.496 to10.884. It cannot recover missing faint
objects or create a mass-complete 1e9 sample at these redshifts. Raw contrast
16.59→12.56pp; paired mass/z mixture9.460→9.647pp. Seven-property raw panels
change the population; only the paired quenching panel holds the mixture
fixed. Parent/cut counts90,952/52,620. These are a wedge, not all-sky results.

## Interpreting the HOD question

The published Abacus BGS mock explicitly adds statistical central/satellite
galaxies, including faint objects and a particle-based extension for unresolved
low-mass hosts. It does not require every satellite to be a resolved subhalo.
Low stellar mass is not low host-halo mass: a faint satellite can inhabit a
massive resolved host. Below the resolved-host limit, the published extension
has negligible satellite occupation; this does not imply satellites are absent
from resolved massive hosts. Brightness selection then excludes sufficiently
faint galaxies. A surviving stellar-mass cut cannot undo this loss.
Source: https://arxiv.org/html/2312.08792v2 (sections2,3.3,5.2).
This establishes the construction, not the central/satellite completeness of
our exact local selected catalogue. Stellar masses/SFR are observed CG15
properties here; the HOD luminosity assignment is not a calibrated model of
the true stellar-mass-dependent environmental quenching relation.

## Planned controlled survey experiments — not executed

Question: how much does each observation step alter the true property relation,
and how much additional error does the environment predictor introduce?
Use one paired mock with persistent galaxy IDs, known real-space positions,
velocities, luminosities/colours, stellar mass/SFR, central/satellite status,
and verified R2/R7 matter truth. Freeze the same spatial label split and
random seeds. Existing TNG snapshot99 supplies a mechanism experiment; its
z=0 properties alone cannot establish evolution at z=.2−.3. A realistic
lightcone requires validated epochs and galaxy/property assignments.

For every stage report TWO quantities: (1) true property relations for the
selected target population; (2) predicted-versus-true relations on those exact
same targets. Also hold a common target subset fixed across stages to separate
losing targets from losing neighbours. Keep fixed environment definitions and
report ranks and physical classes separately, with mass/redshift controls.

1. **Brightness selection.** Assign observer distances and use validated
   apparent magnitudes/K corrections. Apply BGS BRIGHT targeting. Compare
   against random removal matched in redshift and number density, using several
   seeds. The truth change for selected targets measures population selection;
   the extra prediction error measures context loss/bias. Record selection
   fractions by stellar mass and central/satellite status. Freeze properties.
2. **Galaxy motions.** On identical galaxy IDs and brightness-selected targets,
   move positions into observed-redshift coordinates using their velocities.
   Keep physical truth attached to galaxy IDs, rebuild graph inputs, and compare
   against real-space inference. First freeze membership to isolate displacement;
   then separately permit movement across redshift/brightness selection boundaries.
3. **Angular footprint and holes.** Apply the validated mask without fibre
   removal. Compare identical surviving interior and near-boundary queries,
   with/without padded context. Match effective volume when reporting density;
   losing area is not automatically lowering density within retained area.
4. **Fibre assignment / redshift success.** Use the actual targeting and
   assignment rules or a validated emulator, then a separately modelled success
   selection. Compare with random removal matched in angular/redshift completeness.
   Restrict the primary comparison to the same observed target IDs; report
   excluded targets separately. This tests preferential loss of close neighbours.
5. **Interactions and adaptation.** Run brightness × motions × assignment
   combinations under one fixed mask, retaining the sequential practical chain.
   Refit the same model/budget for each condition and also apply one frozen model.
   Hold input transformation conventions explicit; compare frozen training scaling
   and observed-graph refitting as separate arms, without choosing by property
   trend strength. Match graph-scale definitions or add that mismatch as a
   separate factor. Effects depend on ordering; do not sum sequential changes
   into a unique universal causal budget.

Acceptance gates: unique-ID crosswalk, unit/cosmology/coordinate checks, exact
truth attachment under coordinate moves, verified apparent-magnitude selection,
no property leakage, valid common support, paired spatial uncertainty, and
prediction/uncertainty validation against mock truth. Do not use the known
invalid Abacus low-z mapping to fill missing truth. Reserved phases remain
untouched. Run costly graph rebuilds/training only on an explicitly approved
compute allocation. A survey run is not launched by this figure-producing task.

This sequence isolates specified mechanisms within the chosen simulation.
Independent galaxy-formation/occupation models are still needed to test how
well that decomposition applies to DESI. It is not necessary to claim that
the current modest contrast is a bug, or that the historical strong contrast
was wholly false, before doing these tests.
