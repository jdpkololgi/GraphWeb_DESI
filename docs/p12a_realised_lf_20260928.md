# Realised parent luminosity diagnostic — 2026-09-28

Executed workflows/p12a_vac/realised_parent_lf.py on approved CPU allocation
59019148/nid004217. Exposed ph000 only, deterministic stride20, all delivered
parent rows on the existing common mask; no apparent magnitude cut for parent
counts. Stored R_MAG_ABS thresholds, RSD redshift shells, geometric c000 volumes.
No catalogue changes, fitting, new phase exposure or production inference.

## Findings

V1 cumulative density for stored M<-22 relative to recovered wsys LF:

| z | realised / wsys | fraction passing 12<=r<19.5 |
|---|---:|---:|
| .15–.25 | .981 | 1.000 |
| .25–.35 | .983 | 1.000 |
| .35–.45 | 1.003 | 1.000 |
| .45–.55 | 1.029 | .746 |

The realised population therefore broadly reproduces the recovered table in
this well-populated luminosity regime. The alternative nowsys ratios are
.972,.974,.994,1.020; this does not identify the executed HOD branch.

At .45–.55 and stored M<-22:

| | Old | v1 |
|---|---:|---:|
| Delivered density h^3/Mpc^3 | .00013115 | .00011648 |
| Estimated delivered count (sample x20) |214940|190900|
| Fraction passing apparent cut |.34577|.74552|
| Largest sampled apparent r |20.1998|20.0045|

V1 has ~11% fewer delivered galaxies at this stored threshold, yet ~1.92x as
many pass the apparent cut. This is evidence against explaining the selected
count difference solely by a larger intrinsic bright population. Mapping,
colour distributions and luminosity definitions matter strongly. It is not
an isolated K-correction causal experiment: stored M conventions differ and
parent truncation must be respected.

At M<-21 v1 density/table is .750 in the highest shell, where delivered r
reaches20.2; this is consistent with parent flux truncation and must not be
interpreted as failure of the underlying LF. At M<-23 there are only11,18,25,29
sampled v1 objects per shell; apparent ratios1.70,1.39,1.22,1.02 have insufficient
precision for strong bright-tail claims. M<-24 has no sampled objects; undefined
maximum magnitudes are represented as NaN in the Python-readable result.

## Source interpretation correction

Calling the recovered table an executed "input LF" was too strong. Archived
hod_bgs_abacus.py defaults redshift_evolution=False; it loads the target LF
only if evolution is enabled. lookup.py describes the wsys table as predicted
from best-fitting HODs, preserving clustering at z=.2. Agreement therefore
supports consistency with that HOD-predicted LF, not proof of historical
execution settings. Earlier independent published-LF comparison remains a
valid table comparison, but cannot establish an executed LF input or cause.

## Limits and next discriminator

Delivered parent densities are not guaranteed complete LF estimates. Maximum
r below the delivered flux limit is supportive, not proof that every colour
or missing source population is represented. Common mask is empirical; volumes
are geometric, not official effective-volume estimates. One deterministic
sample supplies no formal covariance. No independent Loa luminosities used.

Next decisive test: reconstruct old and v1 reference luminosities under a
single independently calibrated K/E convention, separate colour and K effects
on fixed parent rows with support checks, and compare with Loa luminosity-colour
counts in volume-complete cells. Existing fixed-row mapping swaps and the new
bright-selection fractions prioritise this over a simple global LF normalization
repair. Retain provisional VAC decision; no evidence for changing Loa to fit mocks.
