# Joint colour, luminosity and redshift diagnostic — 2026-09-28

## Result

The surplus in canonical v1 is not solely an overall normalization difference.
A conditional colour mismatch remains at high redshift. On identical fine
z/r cells shared by old parent, v1 and corrected Loa, mean mock-minus-Loa g-r
at 0.45<=z<0.55 is:

| PHOTSYS | Old v0.1 | v1 | Corrected Loa weight retained |
|---|---:|---:|---:|
| S | -0.107 mag | -0.074 mag | 80.0% |
| N | -0.081 mag | -0.119 mag | 61.0% |

Negative means bluer. v1 is not uniformly better across photometric regions.
These are conditional exploratory offsets, not fitted calibration corrections
or independent-realization significance tests. Pair-specific support results
are also retained in RESULTS.json; do not rank versions on different supports.

The prior full-count comparison uses Galactic caps, NOT these PHOTSYS groups:
old/raw divided by corrected Loa is1.087/1.032 at .15-.25 SGC/NGC,
1.021/.981 at .25-.35, .876/.880 at .35-.45, and .553/.589 at .45-.55.
v1 is1.195/1.112,1.194/1.107,1.243/1.193,1.508/1.561 respectively.
Thus the old discrepancy begins before the final bin and steepens in the tail;
v1 has a broader surplus that also grows toward the tail. More parent galaxies
provide selection headroom, but do not establish the correct conditional
population or galaxy-matter relation.

## Method and validation

Deterministic stride20 from canonical exposed ph000 v0.1/v1 and Loa v2.1;
1,013,159 selected rows total. Common NSIDE256 sky, numerical12<=r<19.5,
.15<=z<.55. Loa ZWARN0, DELTACHI2>=25, SPECTYPE GALAXY; positive g/r flux;
weights WEIGHT_ZFAIL/(FRACZ_TILELOCID*FRAC_TLOBS_TILES). No FKP or imaging
fitting weights. PHOTSYS from Loa; recovered DEC32.375 split for parents.
Histograms dz=.01,dr=.1,d(g-r)=.05. Conditional comparison uses cells with
>=20 sample/weighted counts in all three sources and standardizes mock cell
weights to Loa ONLY for this diagnostic. No catalogue or selection changed.

M_r-5logh is reconstructed for all samples with the same recovered v1 colour
inversion, r K-correction and Q=.67 evolution, approximate native flat-LCDM
luminosity distance. This tests a common-mapping hypothesis. It is not an
independent luminosity measurement or volume-corrected LF. The high-z shape
plot retains differences after unit normalization; magnitude and colour are
coupled by the adopted transformation, so do not count them as independent
causal evidence. Exact imaging/quality parity remains open.

Validation: FITS memmap/fitsio equality on100 strided rows per source;
finite/nonnegative histograms; mock counts conserve; colour and M histogram
weights conserve; no out-of-range selected colour/M rows. Stride-expanded Loa
weight is1.000479 times the existing full-scan value. Deterministic sampling
is not random sampling and no naive Poisson significance is claimed.
Initial fitsio strided reader was I/O-bound; stopped step59010842.2 and used
validated memmap reads. Fine-bin run and plotted shared-support summaries pass.

Evidence: evidence/p12a_joint_colour_luminosity_20260928/{RESULTS.json,
HISTOGRAMS.npz,VALIDATION.json}; figure figures/p12a_joint_colour_luminosity_20260928/joint.png.
Scripts: joint_colour_luminosity.py, summarize_joint_colour_luminosity.py,
plot_joint_colour_luminosity.py in workflows/p12a_vac, run in that order.

## Next discriminating work

1. Complete pinned imaging/assignment/veto replay before comparing selected
   populations; preserve independent spectroscopic-success limitations.
2. Repeat conditional comparison after those stages, retaining PHOTSYS split,
   bright/faint magnitude dependence and reported common-cell support.
3. Check luminosity against independent fastspecfit/DR2 LF conventions and
   published Moore et al.2511.01803v2 tables before attributing discrepancies
   to HOD/LF or K/evolution. A common v1 transformation alone cannot decide this.
4. Require clustering and galaxy-matter validation before any replacement VAC.

Operational correction: allocation59010842 was RUNNING but original forward
step.0 had CANCELLED by0 / signal9 after2m43s. No application traceback or
completed preparation was present; exact cancellation cause remains unproven.
It was reused for this bounded joint diagnostic. Allocation state must not be
reported as an active scientific pipeline without checking its step/process.

Preparation handoff: interactive two-row mask smoke passed on59010842 with
finite mask/NOBS values and row conservation. Full preparation previously
reached32-worker configuration but was cancelled before completing input read.
Submit only frozen preparation+validation via regular batch,64CPU/one node/2h;
32 workers bound memory/I/O concurrency on the exclusive node. Source checksums
and LSS revision are checked at launch. No auto retry; no untested ledger or
full replay batch chain. CFS read-only; isolated existing Scratch output,
overwrite disabled. Batch success requires PREPARATION_READY.json, not Slurm
COMPLETED alone. This handoff uses existing user compute authorization.
