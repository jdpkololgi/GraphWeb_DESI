# Historical BRIGHT footprint test — 2026-09-28

## Result

The supplied Ashley/Arnaud thread identifies a real geometry mismatch. Replaying
Ashley's exact historical tile path and BRIGHT & (IN_DESI==1) selection does
**not** explain our common-sky high-z deficit. It removes only 0.330% of SGC and
0.192% of NGC common-mask area. The old-parent high-z deficit persists in the
matched footprint and an eroded interior; v1's surplus also persists.

Historical tiles:
`/global/common/software/desi/perlmutter/desiconda/20230111-2.1.0/code/desimodel/main/data/footprint/desi-tiles.ecsv`.
IN_DESI is int16 (0/1), not bool; 5994 tiles selected. Explicit comparison avoids
integer-mask indexing ambiguity. File hash is pinned in RESULTS.json.

## Method and validation

Script: workflows/p12a_vac/historical_footprint_test.py.
CPU allocation59018519/nid004175, deterministic stride20 scans of exposed ph000
old and v1 raw parents and Loa v2.1 full_HPmapcut. Same .15<=z<.55, numerical
12<=r<19.5 and existing NSIDE256 common mask. Loa ZWARN0, DELTACHI2>=25,
SPECTYPE GALAXY; weights WEIGHT_ZFAIL/(FRACZ_TILELOCID*FRAC_TLOBS_TILES).
No mock quality observables invented; raw parents versus corrected observed
population remains an imperfect observation-model comparison.

Each selected object has exact is_point_in_desi membership in historical tiles.
For all347202 sampled selected old-parent galaxies, geometry agrees exactly
with stored IN_Y5:346364 both inside,838 both outside,zero disagreements.
This verifies the historical restriction in this sample, not every forFA file.

Three domains: existing common mask; its intersection with exact historical
footprint; intersection with an interior mask formed by twice removing pixels
adjacent to uncovered NSIDE256 pixels. Final interior object membership also
requires exact historical membership. Two pixel rings are an operational edge
exclusion, not a fixed angular-distance cut. No catalogues/randoms altered.

Equal-area HEALPix quadrature at NSIDE512 and1024 estimates areas; maximum
relative area difference0.006%. Shell volumes use c000 native approximate
Omega_m=(.02237+.1200+.00064420)/.6736^2, comoving distances in Mpc/h, volume
Omega*(chi_hi^3-chi_lo^3)/3. RESULTS.json provides fine-bin densities and counts.
These are common **geometric** volumes, not newly certified official random-
derived effective volumes including every subpixel imaging/veto selection.
The same denominator is used for mock and Loa, so their ratios do not depend
on area-estimator normalization. Source size/mtime are checked unchanged.

| Domain | SGC area deg2 | NGC area deg2 |
|---|---:|---:|
| Existing common |4176.955|9543.870|
| Historical intersection |4163.159|9525.567|
| Interior |3395.572|7967.415|

## Counts / corrected Loa (same geometric volume)

| Population / domain | .15–.25 S/N | .25–.35 S/N | .35–.45 S/N | .45–.55 S/N |
|---|---|---|---|---|
| Old, original common |1.094/1.032|1.016/.976|.875/.887|.540/.582|
| Old, historical matched |1.093/1.031|1.016/.976|.873/.887|.541/.582|
| Old, interior |1.056/1.010|.984/.960|.836/.864|.532/.582|
| v1, original common |1.190/1.111|1.188/1.098|1.225/1.213|1.535/1.547|
| v1, historical matched |1.190/1.111|1.187/1.098|1.223/1.213|1.534/1.545|
| v1, interior |1.153/1.088|1.152/1.083|1.178/1.190|1.492/1.537|

S/N here means Galactic SGC/NGC, not photometric N/S. Deterministic stride
explains differences from prior full-count estimates; geometry contrasts use
identical sampled objects. No covariance/significance claim from one phase.

## Consequence

Geometry is not a sufficient explanation for our residual; photometry/LF and
population investigation remains necessary. This does not prove a specific LF
cause, exonerate all observation effects, or measure low-k clustering impact.
Preserve the historical footprint in any old mock/random comparison. No grounds
here to remove additional Loa galaxies from the primary provisional VAC.
Formal official-random effective-area normalization and full forFA/mock-random
parity across phases remain separate checks; this task establishes the much
narrower conclusion that the historical restriction does not remove our count
trend on already-common sky. Existing provisional-release decision unchanged.
