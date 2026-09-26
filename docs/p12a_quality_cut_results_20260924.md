# P12-A DESI quality-cut census, 2026-09-24

This read-only census scanned every row of both catalogues under CPU allocation
58824338. Evidence: `p12a_quality_cut_audit_20260924.json`; implementation:
`workflows/catalog/p12a_quality_cut_audit.py`. No catalogue was modified.

The historical GraphWeb builder requires ZWARN==0, DELTACHI2>=25,
SPECTYPE==GALAXY and BGS BRIGHT membership. Installed LSS common_tools.goodz_infull
for BGS requires ZWARN==0 and DELTACHI2>40; it does not itself impose GALAXY.
The LSS helper alone is not an authority for replacing an already validated
scientific sample. A difference is acceptable if the selected population and
response are modeled consistently. The earlier mismatch warning was too broad:
it did not demonstrate that the existing DESI catalogue violated its own cuts.

| Census, 0.15 <= z < 0.55 | Existing GraphWeb >=25 GALAXY catalogue | Loa full HPmapcut candidate |
|---|---:|---:|
| All rows in this redshift range | 6,466,730 | 5,514,226 |
| ZWARN==0 | 6,466,730 | 5,504,897 |
| Historical >=25 and GALAXY rule | 6,466,730 | 5,436,413 |
| Installed LSS >40 helper rule | 6,449,978 | 5,481,423 |
| Historical rule only | 16,752 | 13,293 |
| LSS helper rule only | 0 | 58,303 |

All9,166,391 rows of the existing GraphWeb catalogue satisfy its advertised
ZWARN, DELTACHI2 and SPECTYPE cuts. Tightening its threshold removes0.259% in
the VAC range, rising to1.019% in the highest shell. The Loa full table's two
rules are not nested: LSS admits high-confidence non-GALAXY spectra while
rejecting the25-to40 band. These are separate effects. The different source
footprints/parent populations also mean the two table columns are not an
identity-matched comparison.

P12 ph002-005 mock observed-success selection is finite positive Z_not4clus
and ZWARN==0. It does not impose measured DELTACHI2 or SPECTYPE: mocks supply
true simulated galaxies rather than noisy spectral classifications. Consistency
therefore requires modeling the real successful-galaxy population/response;
it cannot be established by copying spectroscopy flag cuts onto absent mock
columns. Neither changing to40 nor retaining25 automatically requires training
again. Preserve both possible real cuts as explicit diagnostics; freeze the
chosen Loa observed sample and its count/support response together before the
real-data adapter is qualified. Do not substitute the historical zall catalogue
for the Loa full catalogue without a footprint/identity crosswalk.
