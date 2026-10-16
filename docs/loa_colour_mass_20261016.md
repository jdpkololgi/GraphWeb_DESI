# Observed colour–mass diagrams for the saved Loa selections

Recreated 2026-10-16 directly from the frozen selected FITS products, without
rerunning selection, completeness estimation, inference or population matching.

| Redshift | Minimum log10(Mstar/Msun) | Selected | Red | Blue | Outside colour cuts |
|---|---:|---:|---:|---:|---:|
| 0.15–0.20 | 10.60 | 534,070 | 518,962 | 14,909 | 199 |
| 0.20–0.22 | 10.69 | 232,024 | 229,932 | 1,922 | 170 |
| 0.20–0.30 | 10.96 | 611,197 | 532,770 | 375 | 78,052 |

All intervals include the lower redshift endpoint and exclude the upper.
Red is strictly 0.75 < g−r < 1.5; blue is strictly g−r < 0.75.
Colours are observed, Milky-Way corrected, without K correction. Colour
membership is not a quenching label and changes with redshift.

The three panels share axes, bin edges (0.0125 dex by 1/60 mag) and logarithmic
count normalization. Grey shading marks masses below the fixed selection;
vertical dashed lines mark completeness cuts. Horizontal lines mark g−r=0.75
and 1.5. The distributions include selected galaxies outside the colour cuts.
The common display bounds omit 1, 1 and 7 finite objects respectively; these
counts are annotated, and no selected object has nonfinite plotted values.
These are raw counts, not volume-normalized densities or matched populations.

[Combined PNG](figures/loa_colour_mass_20261016/colour_mass_comparison.png) ·
[Combined PDF](figures/loa_colour_mass_20261016/colour_mass_comparison.pdf).
Individual PNG/PDF figures are available in the same directory with suffixes
`z015020`, `z020022`, and `z020030`.

Source: `workflows/presentation/loa_mass_complete_20261016/colour_mass_diagrams.py`.
Evidence: `docs/evidence/loa_colour_mass_20261016/MANIFEST.json` and
`histograms.npz` retain input hashes, code hash, exact bins, counts and outputs.
The approved CPU allocation was 59858009 on nid200041; render exited zero.
Allocation release was verified by an empty squeue result with no stderr.
Checks passed for unique IDs, redshift and mass membership, exact saved
red/blue counts, histogram sums and independent flat-bin reconstruction.
The combined PNG was visually inspected; repository copies are byte-identical.

The repository graph refresh was not run: graphify was unavailable, and automatic
approval review separately rejected `graphify update .` for unspecified mutation
and network effects. Figure generation and scientific validation were unaffected.
