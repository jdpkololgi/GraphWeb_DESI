# Completeness-corrected Loa and pinned v1 forward run — 2026-09-28

## Completed comparison

Loa full_HPmapcut v2.1, common NSIDE256 sky, ZWARN=0, DELTACHI2>=25,
SPECTYPE=GALAXY, .15<=Z_not4clus<.55. 5,397,649 selected rows; zero invalid
completeness factors. Use 1/(FRACZ_TILELOCID*FRAC_TLOBS_TILES), with delivered
WEIGHT_ZFAIL as a separate sensitivity. No FKP, imaging-density or n(z) fit.
This corrects assignment for the selected spectroscopic population; it is not
proof of complete underlying imaging selection or absolute redshift success.
The delivered failure weights model relative trends, not necessarily all failures.

Under common12<=r<19.5, rawv1 / correctedLoa including WEIGHT_ZFAIL:

| z | SGC | NGC |
|---|---:|---:|
| .15–.25 | 1.195 | 1.112 |
| .25–.35 | 1.194 | 1.107 |
| .35–.45 | 1.243 | 1.193 |
| .45–.55 | 1.508 | 1.561 |

Using all historically targeted Loa objects gives high-shell1.508/1.479,
but rawv1 still has a uniform19.5 cut in that particular bracket, so North
selection is not identical. Raw parents have common pixel support, not exact
subpixel imaging/veto masks. These findings motivate the actual forward run;
they are not final matched-mask n(z) estimates. Corrected Loa weights cannot
be substituted for an unweighted observed point process in environment training.

Evidence: `evidence/p12a_loa_completeness_20260928/RESULTS.json` and code
`workflows/p12a_vac/loa_completeness_counts.py`. Uniform observed counts replay
the earlier full Loa magnitude histogram exactly. Actual installed LSS code at
commit d942b990860e016e7558c740966fea3e5ce7e7f3 implements this product of
completeness factors in py/LSS/main/cattools.py (weighttileloc,compmd=dat).

## Broader search: useful implementation evidence

- Smith et al. https://arxiv.org/abs/1809.07355: BGS fibre assignment depends on
  projected density. Historical retention cannot be assumed portable to denser v1.
- https://arxiv.org/html/2411.12020 section11.2: DR1 altMTL reproduces assignment
  history/hardware. Its mock redshift failures are imposed by random rejection
  calibrated to sample good-redshift counts, not synthetic DELTACHI2. This is
  useful implementation precedent, NOT evidence for Loa high-tail agreement or
  authorization to tune our mock counts to data.
- Direct GitHub tree inspection of desihub/LSS located DA2 BRIGHT preparation,
  initializer, altMTL loop and Loa catalogue scripts. Pinned source snapshots
  retained in evidence/p12a_v1_forward_20260928/source.

## Forward-run contract and execution

One exposed ph000 canonicalv1; scratch root
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/v1_forward_ph000_20260928/`.
No writes to CFS, no other phases. Upstream LSS source pinned d942b990.
Installed compute-node LSS hash agrees; desitarget5.5.0.post6118,
fiberassign6.0.3.post3871 under desiconda20260910-3.1.0.

Preparation adapter: v1 input schema (no fabricated IN_Y5), seed20260928,
BRIGHT North/South limits19.54/19.5 at recovered DEC32.375split,
standard upstream FAINT competitor prescription (69.5%sampling,20%promotion),
Y3 tiles and imaging-mask application, no BGS bright density downsampling.
Globals use loa-v1 rather than upstream iron. Upstream ZWARN0 is explicitly
an initial ledger placeholder; it does not simulate spectral quality.

First attempt exposed upstream unguarded multiprocessing Pool incompatibility
with Python3.14 forkserver. Single-phase adapter executes serially instead.
Initializer adaptation limits workers to8 and pins Linux fork start method.
This is an implementation repair, not a population-model adjustment.

Full forward outcome is pending. Preparation is not fibre-assignment completion,
and no success is claimed for absent quality observables or unexecuted stages.

Operational follow-up: stopped step59009604.5 before the upstream masking
Pool(default128workers) would hit the same forkserver recursion. Current
adapter pins fork and8workers in both masking calls; rerun log prepare_fork.log.
Initial/second attempt logs retained; no scientific setting changed by this fix.

Further source lead: https://astro.dur.ac.uk/~cole/BGS/LFs.html supplies Moore
et al.(2026) North/South g,r,z,w1 luminosity functions (Vmax and SWML), including
red/blue subsamples. These are an independent published population-model check,
not evidence that the November2024 canonical file used the published tables.
No replacement LF or calibration to Loa was applied.

Comparison figure: docs/figures/p12a_loa_completeness_20260928/counts.png/pdf,
created by plot_loa_completeness.py using canonical plot style and inspected.
Uniform observed histogram replay against previous evidence passes exactly.
