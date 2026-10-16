# P12-A CG5 versus CG15 check

All seven property-environment trend directions persist with CG5. The fits use
the same P12-A probabilities and the same jointly valid galaxies. The larger
legacy-GAT versus P12-A contrast gap is not explained by choosing CG15 here.

There are 102,621 jointly valid galaxies in the existing z=.20-.30 wedge.
Common mass/redshift support retains 75,321 for non-mass properties.

| Cluster-minus-void contrast | CG15 | CG5 |
|---|---:|---:|
| Probability-weighted quenched fraction | +8.970 pp | +9.573 pp |
| Probability-weighted log sSFR | -0.22736 dex | -0.21773 dex |
| Hard-label quenched fraction | +11.764 pp | +12.291 pp |
| Hard-label log sSFR | -0.29414 dex | -0.28175 dex |

CG5 shifts overall SFR/sSFR levels upward but leaves similar environmental
contrasts. Using CG5's own masses for mass/redshift controls gives +8.894 pp
and -0.20320 dex (probability weights; 78,317 retained rows). This sensitivity
changes support and reference weights and is not a strictly paired contrast.

Primary controls: .1-dex CG15 mass/.01-redshift cells, same reference mixture
for both fits and both membership methods, at least 20 effective objects in
every class. Stellar mass controls redshift only. Quenched: log sSFR<-11.
Paired 64-draw sky-block resampling gives CG5-minus-CG15 change +0.603 pp
[+0.437,+0.793] and +0.00963 dex [+0.00771,+0.01200] (16-84% intervals).
These small changes are measurable; they do not explain away the trends.

Exact IDs, redshifts and probabilities match the existing cache. All ten
CG15/CG5 properties match fresh source reads for 64 sampled galaxies. No new
inference or full source scan. These fits share observations and assumptions;
this is robustness to the fit variant, not independent truth validation.

Source: GraphWeb_DESI/workflows/presentation/p12a_cg5_check_20261014.py.
Machine-readable results: RESULTS.json. Figure: cg5_vs_cg15.png/pdf.
