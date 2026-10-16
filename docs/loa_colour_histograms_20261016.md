# Observed g−r distributions by inferred environment

The original three CIGALE mass selections are retained: 0.15≤z<0.20 at
log10(M/Msun)≥10.60, 0.20≤z<0.22 at ≥10.69, and 0.20≤z<0.30 at ≥10.96.
Colours are the same Milky-Way-corrected observed Legacy Survey g−r as the
preceding colour–mass diagrams. No FastSpecFit masses or rest colours enter.

[Highest-probability classes](figures/loa_colour_histograms_20261016/observed_colour_hard.png)
are the primary view; a [probability-weighted version](figures/loa_colour_histograms_20261016/observed_colour_soft.png)
retains uncertain class assignments. Both have PDF counterparts. Each figure
shows all four environments overlaid, with a linear density view above and
the same histogram on a logarithmic density scale below. Zero bins are absent
from the logarithmic view, not replaced with pseudocounts.

The script `workflows/presentation/loa_rest_colour_20261016/colour_histograms.py`
marginalises the previously validated `diagram_bins.npz` over mass and merges
adjacent colour bins to width 1/30 mag. Each environment is independently
normalised to unit integral over saved display bounds (10.3≤logmass≤12.5,
−0.5≤g−r≤2.5). Those bounds omit 1, 1 and 7 galaxies, respectively, from the
original selections; corresponding probability-weight omissions are recorded.
There is no further red/blue colour cut. The line at 0.75 is a reference only.

All 24 marginal distributions reconcile to the archived displayed counts;
independent column-pair sums agree and all density integrals equal one.
CSV bins, source/script/output hashes, totals and omissions are in
`docs/evidence/loa_colour_histograms_20261016/`. Both PNGs were visually checked.
Only the 462 KB aggregate was read; no catalogue processing or Slurm job was
needed for this brief single-thread rendering.

The plots show a dominant red peak with a blue-side shoulder/tail, rather than
a clearly separated second peak in these mass-limited observed-colour samples.
Clusters have a more concentrated red peak and voids a larger relative blue
shoulder. These are descriptive raw distributions: different mass/redshift
mixtures and selection can contribute, no matching or volume weighting has
been applied, and the higher-redshift windows overlap. No causal environment
claim or formal bimodality test is implied.

The repository's `graphify` command is unavailable; targeted source inspection
was used instead.
