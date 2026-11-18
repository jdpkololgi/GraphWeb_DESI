# Loa density structure and tidal sensitivity audit

## Scope and status

Existing NGC expanded U-Net seed17, four shared-state density draws; no new
inference, training, mock truth access, VAC replacement or website publication.
Source: `workflows/presentation/loa_structure_audit_20261028.py`.
Final outputs: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_structure_audit_20261028_v3/`.
`RESULTS.json` binds source/operator hashes; `diagnostic.png` and
`tidal_sensitivity_layers.npz` provide provisional local layers, not a validated
full-survey tidal release. The existing 80-second video is unchanged.

## What the display mask means

The video renderer uses `(support >= 0.5) & (coverage > 0)` and rejects
interpolation unless all eight corners pass. This is a display convention,
not a calibrated inference-quality threshold or a physical void criterion.
Native support is the average of eight binary finer-grid footprint cells,
so 0.5 means at least four of eight supported finer cells. Coverage separately
counts prediction windows. Support is not the number of observed galaxies.

The region-selection rule uses a 65-cell uniform-filter neighbourhood with
more than 98% passing cells, then chooses the candidate nearest the supported
geometry centroid, independently of density. It does not demand a hole-free
region. Neither the 0.5 cutoff nor the 98% threshold was calibrated against
scientific accuracy; both are visualization heuristics. Their operational
purpose is to avoid showing weakly supported/edge extrapolations on equal terms
with observed-footprint predictions. They do not establish that excluded
predictions are invalid, nor that included predictions are trustworthy.

Measured central 32^3 core: all 32,768 cells have prediction coverage and finite
density in all four draws. Exactly 680 (2.0752%) have support below 0.5;
241 have zero support. Retained targets: 32,088. These are mask holes, not failed
inference. Separating genuine angular footprint holes from rasterization effects
requires a further mask audit; it was not established here.

## Structure and resolution

Frozen target source `e2e_coupled_operators.smoothed_spectrum` applies
exp(-k^2 * 7^2 / 2) to native density before target extraction. R7 is therefore
a training-target property, not added by the renderer. Fine target cells are
6.766 Mpc/h; the video overview additionally averages 2x2x2 cells to 13.532 Mpc/h.
Trilinear interpolation and line-of-sight alpha compositing further affect
appearance, but their individual contributions have not been quantified here.

In the fixed 88^3-cell region, on identical retained spatial support, native
density standard deviations across the four draws are 0.4595–0.4691; after
2x pooling they are 0.4185–0.4276 (about 9% lower). No new smoothing is applied.
The native slices show peaks and extended connections. There are 1,045–1,102
3x3x3-grid local maxima above twice the mean density per draw, excluding a
two-cell mask boundary. These are grid peaks, not a halo or cluster catalogue.
This supports the presence of structure in the generated smoothed fields,
not recovery of true DESI filaments, resolved halo interiors, or thin filaments.
Displayed pooled slices average two depth cells; they are not identical-thickness
versions of the native one-cell slice. Masked pixels are explicitly blank.

## Tidal checks and limits

Used the existing cuboid-integrated free-space potential-Hessian operator,
without extra smoothing, evaluated in a fixed central 216.512 Mpc/h cube.
Source cubes span 649.536, 974.304 and 1299.072 Mpc/h. Unmodelled exterior
is assigned mean density, not fitted or sampled. Two source assumptions are
compared: all prediction-covered cells, and fractional footprint occupancy.
The latter changes the source field and hence the trace near partial support;
it is a sensitivity test, not a correction or an alternative posterior draw.

Direct summation versus FFT: max absolute error 6.66e-16; changed FFT padding:
2.78e-16; isotropic cube self term: 1.11e-16. Trace identity errors on retained
targets are at most 1.18e-14. Classes count ordered eigenvalues above 0.2:
0 void, 1 wall, 2 filament, 3 knot. Unretained class cells carry sentinel 255.

Increasing source width from 974 to 1299 Mpc/h changes 0.688% of covered-source
labels and 0.544% of fractional-support labels across the four draws and retained
cells. Changing the support assumption at the largest width changes 2.619%.
These small aggregate changes do not certify each cell or unknown external tides;
the context increments are not monotonically decreasing, so no convergence claim
is made. Local filament/knot labels are present but remain closure-dependent.
No calibrated real-DESI uncertainty or mock-truth accuracy follows.

## Execution and retained failures

User approved one CPU node, 32 CPUs, desi/interactive, Scratch, at most 60 minutes,
with cleanup. Allocation 60418676 ran on nid200008. First attempt stopped on an
overstrict all-supported-core assertion before large convolutions; v1 retained.
v2 explicitly masked unsupported targets. v3 additionally matches spatial masks
for pooled/native dispersion and corrects the pooled-display mask; v2 retained.
The original fields and release are unchanged. Scientific validation beyond these
numerical/sensitivity checks remains open; do not promote these layers as validated.

Allocation released after 6m21s; accounting reports COMPLETED 0:0. Final diagnostic
step 4 completed 0:0. Required local Graphify refresh (step 5) failed with the
existing fcntl.flock OSError 524; no bypass or global-index update was performed.
