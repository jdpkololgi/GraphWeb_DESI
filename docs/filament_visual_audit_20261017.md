# Projected strings versus saved filament probabilities

The observed mismatch is real in the plotted catalogue values, but its cause is unresolved. Exact frozen VAC row/ID/coordinate/support/quality/P(filament) checks passed for all4,459 original galaxies. Of734 targets in the upper20% of projected neighbourhood elongation,482 have P(filament)<0.2:297 wall-argmax and185 void-argmax. This is not a measured model error rate. The statistic is based on the same sparse observed points and is not tidal ground truth.

The canonical T-web label counts eigenvalues above a threshold at a specified smoothing scale. It need not match a fine galaxy string in projection: a slice through a sheet can appear linear, off-slab matter affects the 3D tidal environment, and redshift-space positions include peculiar velocities. These are candidate explanations, not findings established for these galaxies. Different web methods classify different aspects of structure; see [Forero-Romero et al.2009](https://arxiv.org/abs/0809.4135) and [Libeskind et al.2018](https://academic.oup.com/mnras/article/473/1/1195/4062204). Model transfer or calibration failure remains possible.

## Nested-slab screen

Use fixed central907 targets with |height|<2.5Mpc, original x/y edge margin, radius10.34585Mpc, >=8neighbours and the original fixed elongation threshold0.584773. Probabilities never change.

| Full depth, Mpc | Shape-valid targets | Line-like targets | Low-Pfil among line-like |
|---|---:|---:|---:|
|5|454|176|39|
|10|699|214|81|
|20|907|157|102|

Restricting to the454 targets valid at all three depths gives line-like counts176/105/31 and low-Pfil counts39/15/9. Thus93 of102 low-Pfil line-like targets in the20Mpc central sample fail the thinnest-slab neighbour-count requirement. The apparent lower mismatch fraction in thin slices is strongly confounded by eligibility and changing line-like membership. This screen cannot establish or exclude a slab-selection explanation.

Zero complete radius10.34585Mpc spheres fit inside the original20Mpc slab. A decisive test needs a padded off-slab parent sample, common-target thickness/orientation controls, count-matched thinning, complete3D neighbourhoods, and matched held-out mocks with tidal truth. Rotating the already truncated rows would not supply a valid viewing-angle control. No retraining decision follows from this screen.

## Presentation assets

- four_environment_probabilities: identical4,459 rows and common0–1 scales show where probability goes.
- filament_probability_highlight: fixed geometric central200×200Mpc crop,1,539 galaxies retained;626 with Pfil>=0.5 highlighted. This displays filament-probable galaxies, not an extracted filament network. No favorable-region search or retuning.
- local_environment_illustrations: hand-drawn conceptual void/wall geometry, actual40.8/59.2% from512 saved P12-A eigenvalue draws for ph006 parent2717758, and the same242 actual observed-success neighbours. Orientation and local shape are illustrative, not inferred. Scope follows the user's explicit per-galaxy framing.

Package: /global/homes/d/dkololgi/fmsu_tokyo_20261010/filament_and_local_geometry_20261017/.
Scripts: workflows/presentation/desi_filament_visual_audit_20261017.py, desi_slab_sensitivity_20261017.py; Illustris/workflows/presentation/tokyo_local_geometry_20261017.py.
Copy/paste investigation and conditional-retraining prompt: docs/filament_investigation_prompt_20261017.md.
No catalogue mutation, new inference, training or Slurm action. Full catalogue checksum inherited; extracted-row integrity checked afresh. Results are exploratory and real-DESI calibration remains unverified.
