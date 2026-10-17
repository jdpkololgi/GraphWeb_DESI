# Why a tidal wall can have string-like galaxy tracers

A tidal wall is a statement about the local smoothed gravitational field. In
our ordered convention lambda1 <= lambda2 <= lambda3, a wall has exactly one
eigenvalue greater than 0.2: lambda2 <= 0.2 < lambda3. A filament has two.
This is a thresholded proxy for the number of strong compressive directions,
not a direct count of completed collapse axes or a fit to the galaxy positions.
The physical distinction is described by
[Forero-Romero et al. (2009)](https://academic.oup.com/mnras/article/396/3/1815/1747417).

The position covariance used in our audit measures the spread of a finite
sample of observed galaxies around its centroid. It is a different operator
from T_ij = partial_i partial_j Phi_R. Schematically in Fourier space, the latter
weights density modes by (k_i k_j/k^2) W_R(k), up to the adopted gravitational
normalization. It therefore includes the smoothed matter distribution and its
gravity, not only the galaxy positions inside one top-hat ball. Numerically
matching the ball radius to the smoothing scale does not equate these operators.

Several effects can make the tracers look like a string:

1. **Structure exists on multiple scales.** A narrow chain of haloes or a smaller
   ridge can lie within a broader environment that is classified as a wall at
   7 Mpc/h. Bright selected galaxies need not populate a sheet uniformly. The
   importance of scale and the different physical fields used by classifiers
   is discussed in [Libeskind et al. (2018)](https://academic.oup.com/mnras/article/473/1/1195/4062204).
2. **Sparse sampling affects covariance.** The examples have only 13–17 neighbours
   in the main sphere. A few groups can dominate its major axis. The bright,
   successful-redshift tracer sample is not a uniform sampling of all matter.
3. **Projection and redshift-space positions matter.** A sheet can be seen close
   to edge-on or intersected by a narrow slice. Peculiar velocities can stretch
   galaxy separations along the line of sight, while the truth labels remain
   attached to native host environments. These effects are plausible physical
   contributors, not separately measured causal explanations for these cases.
4. **The label has a threshold.** Crossing lambda2=0.2 changes wall to filament
   even if the continuous gravitational field and appearance change only a little.
   The displayed filament error is close to this threshold.

The earlier full-sphere test established that projection alone is insufficient:
true-wall tracer neighbourhoods often remain line-dominated in 3D. It did not
measure which combination of sampling, hierarchical structure and velocities
causes that outcome. A dense-matter/real-space comparison with controlled tracer
thinning would be needed to attribute those contributions. None is claimed here.

## Reproducible examples

All cases are from the already-exposed ph006 development phase at 0.20<=z<0.30,
with the frozen projected-string cut, Pfil<0.2, wall argmax, boundary clearance
and at least eight main-radius neighbours. Within each case category, choose the
median parent ID; no visual reselection. Candidate pools contain 522 true-wall
line cases, 118 true-wall sheet cases and 28 true-filament errors. These three
examples are illustrative and have no population prevalence interpretation.

| Case | TARGETID | True ordered eigenvalues | P(void,wall,filament,cluster) | Full-R L/P/S |
|---|---:|---|---|---|
| True wall, tracer line | 6032112 | -0.1215, 0.1057, 0.3121 | .1328, .6641, .1992, .0039 | .6637/.2723/.0640 |
| True wall, tracer sheet | 3743465 | -0.2624, 0.1285, 0.3893 | .0176, .9180, .0586, .0059 | .0424/.6441/.3135 |
| True filament, wall argmax | 1493453 | -0.2317, 0.2177, 0.2791 | .1523, .6738, .1680, .0059 | .7620/.1887/.0492 |

First case: the 17-neighbour full sphere is line-dominated, but only lambda3 is
above threshold. Many of its nearby tracers also carry true-wall labels. P12-A
correctly ranks wall highest while retaining appreciable filament uncertainty.
Second: the 14-neighbour full sphere is sheet-dominated, although the projected
slab passes the same string cut; P(wall)=0.918. Third: the 13-neighbour full sphere
is line-dominated and the centre is truly a filament; lambda2 exceeds the
threshold by only 0.0177, while P12-A favours wall. It is a real argmax error,
not evidence by itself that the probabilistic model is miscalibrated.

Each PNG/PDF has galaxy distribution, native tidal truth and P12-A argmax in
matching columns. Top row: exact 20 Mpc target-centred projected slab within
R=7/.6766 Mpc. Middle row: orthogonal view of the full sphere of radius R.
Columns within a row contain identical galaxy positions; the two rows use
different, explicitly labelled neighbourhood cuts. Bottom panels show central
true eigenvalues, original posterior probabilities and covariance statistics.
Colours and marker shapes distinguish the classes; prediction opacity increases
with maximum probability. Stars mark the central galaxy.

Native host truth is plotted at the galaxies' observed coordinates. There is no
continuous matter-density or tidal-field image in these figures. A coloured
point is a per-galaxy label, not a reconstructed spatial volume.

## Data, prediction and verification

Extraction retains 133 unique neighbours in 1.5R padded spheres. All 133 have
available truth and predictions. For neighbours, 512 new draws are generated
with CPU reconstruction of the frozen FMPE checkpoint from the saved ph006
encoder summaries and response context. No training, encoder rerun or sealed
phase access occurs. Seven rows with original evaluation probabilities retain
those probabilities, including all three centres; the other 126 use the new
draw probabilities. Fixed torch seed 171018. These are per-galaxy posteriors,
not coherent joint spatial draws.

Dataset/checkpoint/summary/parent/index/truth/cache hashes are checked. Thirteen
neighbour contexts that also occur in the saved selection sample replay to
1e-6. Full-ball counts and covariance at R and 1.5R reproduce the previous audit;
the displayed projected-slab anisotropies independently reproduce to 1e-10.
Central truth and original probabilities match saved arrays. All final figures
visually inspected; validation and file hashes saved in FIGURE_VALIDATION.json.

CPU59904866, nid200014, desi/interactive, 64 logical CPUs, approved 30-minute
limit and Scratch license. Extraction/inference/initial plotting step completed
0:0 in35s, final plot step0:0 in3s. Allocation released after3m33s; unrelated jobs
untouched.

Sources: `workflows/p12a_vac/filament_geometry/{mock_examples,plot_mock_examples}.py`.
Arrays/manifests: `docs/evidence/filament_mock_examples_20261017/`.
PNG/PDF: `docs/figures/filament_mock_examples_20261017/`.
