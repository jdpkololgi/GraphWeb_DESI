# Loa selection sensitivity of the controlled TNG regressors

The historical flux-limit completeness method modestly strengthens some
selected-sample trends, but does not establish recovery of a missing
environment association. On an identical parent/cut mass-redshift reference
mixture, the full-range 90% cut barely changes the contrasts.

This completed calculation uses the four **local-feature MLP regressors**,
not the old GAT. The actual-GAT experiment is separate.

The declared windows are 0.20–0.30, 0.20–0.25 and 0.25–0.30. These reuse
available predictions; they do not reproduce the historical 0.01–0.06
selection. All five models use identical galaxies within each panel.
Within each dz=0.01 slice, use the faintest 20% relative to the BGS BRIGHT
targeting limit and the 90th percentile of logM+0.4(r-r_limit). The sensitivity
cut adds 0.2 dex. The parent for estimating the limit has no stellar-mass
truncation; outcome samples retain the established valid-CG15 and geometry
requirements and 9<=logM<12. This is a flux-limit completeness diagnostic,
not correction for fiber assignment, surface brightness or SED errors.

Environment trace quartiles are frozen on the uncut bright parent, separately
for each model and dz=0.01 slice. They are not reassigned after a cut. The
control cells are 0.1 dex in stellar mass and 0.01 in redshift; all models
share the reference mixture within a panel. Different panels need not have
the same reference mixture. The mass-property panel controls redshift only.

Full-range mass/z-controlled Q4-minus-Q1 quenched contrasts, percentage points:

| Model | Uncut | 90% cut | 90% + 0.2 dex |
|---|---:|---:|---:|
| Abacus P12-A | 6.38 | 8.34 | 8.24 |
| TNG R2 normal | 4.98 | 6.15 | 6.45 |
| TNG R2 thinned | 8.36 | 10.33 | 10.61 |
| TNG R7 normal | 5.71 | 7.25 | 7.07 |
| TNG R7 thinned | 6.29 | 8.09 | 8.27 |

The corresponding parent sample has 90,952 galaxies; the cuts retain
52,620 and 29,447. Controlled common-cell counts are separately recorded
and displayed in the figures. The one-row difference from the preceding
90,953-row analysis is the explicit apparent-magnitude limit.

Holding the parent and selected samples to **the same** full-range reference
mixture changes the 90% comparison to:

| Model | Parent | Cut |
|---|---:|---:|
| Abacus | 8.669 | 8.552 |
| TNG R2 normal | 6.075 | 6.356 |
| TNG R2 thinned | 10.559 | 10.421 |
| TNG R7 normal | 7.590 | 7.624 |
| TNG R7 thinned | 8.331 | 8.329 |

Thus much of the apparent increase in the first table is associated with
changing population/support. Narrow redshift windows and the stronger cut
show some variation; all declared results are retained, with no choice of a
preferred cut based on observational trend strength. These rank contrasts
must not be equated with the historical cluster-minus-void contrasts.

690 independent direct-mean reconstructions pass, along with exact TARGETID
parity, nested-cut checks and input hashes. Intervals use 128 paired sky-block
resamples (16–84%); they exclude fit uncertainty and transfer/systematic error.
The 21-page PDF contains all seven properties, raw/controlled rank curves
and raw physical-class curves. Empty classes remain gaps.

Code: `workflows/presentation/loa_selection_sensitivity_20261015.py`.
Artifacts: `~/fmsu_tokyo_20261010/loa_selection_sensitivity_20261015/`:
`selection_overlays.pdf`, six summary PNGs, and `RESULTS.json` with full curves,
support, completeness limits and paired-cut comparisons. Execution used the
explicitly authorized one-thread login fallback after Slurm failed both
inside and outside the sandbox. No scheduler job or production VAC changed.
