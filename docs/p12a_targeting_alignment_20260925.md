# Exact selection alignment audit — 2026-09-25

**Exact alignment has not been established.** Previous component ablations
were insufficient to claim it. The user requested causal selection/modelling
investigation, not an imposed correction. This audit supersedes any broader
interpretation of “cuts checked” in earlier reports.

## Evidence and limitations

New full-table Loa/ph006 census under allocation58871893, cosmic_env. Common
sky is the same intersection as the five-phase count plot. Source tables:
Loa-v1/v2.1 BGS_BRIGHT_full_HPmapcut and ph006 observed-with-truth, whose upstream
source is AbacusSummitBGS_v2/altmtl6/kibo-v1/mock6. These different processing
labels do not prove different effective selection, but require a crosswalk.
Source formulas and hashes, available columns, FITS headers and counts are in
[evidence](evidence/p12a_targeting_alignment_20260925.json).
Implementation: workflows/p12a_vac/targeting_alignment_audit.py and
summarize_targeting_alignment.py.

| Component | DESI | Mock | Alignment status |
|---|---|---|---|
| Scientific redshift range | .15≤Z_not4clus<.55 | Same range on successful RSD redshift | Verified downstream; context convention separately frozen |
| BRIGHT bit | Inherited real targeting | Preparation assigns bit2 from R_MAG_APP threshold | Membership passes; bit semantics alone insufficient |
| Apparent magnitude | Extinction-corrected Legacy Survey flux; nominal19.5 south/19.54 north PHOTSYS | R_MAG_APP<19.5 | Numerical difference confirmed; passband/photometry equivalence unproven |
| Colour/fibre selection | Broad g−r/r−z window, fibre magnitude and fibre-total flux limits | No matching fibre-flux or multi-band forward model in inspected forFA product | Not equivalent by construction |
| Gaia star rejection / SGA recovery | Gaia G−r or missing Gaia, with SGA recovery branch | Mock true galaxies, no Gaia or REF_CAT inputs | Full counterpart not established |
| Imaging | NOBS, positive flux inverse variance, mask bits, HP maps | Later angular masks; not noisy galaxy photometry | Some row tests verified; exact production map/version crosswalk open |
| Spectroscopic success | ZWARN0, DELTACHI2≥25, GALAXY | ZWARN0 on noiseless true galaxies; no DELTACHI2/SPECTYPE | Effective success model not demonstrated equivalent |
| Assignment/usable locations | Loa full catalogue | altMTL/kibo full catalogue | Assigned/usable bits checked previously; exact historical recipe not pinned |
| Parent population | Real luminosities/colours | Model luminosities/colours from fixed-snapshot mock | New conditional colour discrepancy; model lineage incomplete |

The installed desitarget/main formulas were inspected and hashed, but are not
claimed to be the exact historical targeting revision. Loa full table lacks
REF_CAT, Gaia G and targeting-release provenance. The FITS headers do not supply
that revision; the inspected v2.1 log.txt is empty. Raw CutSky header gives area
metadata but not generator commit, luminosity-function parameters or passband
recipe. Local preparation scripts describe the path but are not an immutable
receipt for the official ph006 generation. These are specific provenance gaps,
not evidence that the historical versions necessarily differ.

## What the additional row checks show

All5,397,649 successful Loa galaxies on the common sky have the BRIGHT bit.
Almost all also pass the available nominal magnitude, broad colour, fibre
magnitude, fibre-total flux, positive inverse variance, imaging bits1/13 and
NOBS checks; all pass TSNR2_BGS>1000 and have maskbit11 clear. The conjunction
fails at most.002% in any cap/shell (rounding to.001%). Individual failures
cannot be called targeting errors without SGA recovery and original photometry.

The uniformr<19.5 diagnostic removes5.317% of NGC high-shell galaxies and
.002% of SGC. DELTACHI2>40 removes1.364%/1.321% respectively. Applying these
jointly with all available targeting components still leaves DESI/ph006
high-shell count ratios **1.9155SGC / 1.7806NGC**. Thus these measurable
extra/looser real-data cuts do not explain most of the excess.

This is a retained-sample audit, not a reconstruction of rejected targets.
It cannot determine selection completeness from survivors alone. Missing
fibre/colour cuts in a mock would, by themselves at fixed underlying population,
retain more mocks, not explain a deficit. Population-dependent photometry or
selection mappings remain plausible causes. SGC is PHOTSYS=S here, while NGC
mixes N/S: galactic cap and photometric system must not be conflated.

## New upstream population diagnostic

Compared Loa de-extincted g−r with mock G_R_OBS. Both restricted to r<19.5;
within each cap, mock colour histograms are descriptively reweighted to Loa's
counts in Δz=.01 bins before combining into .1-wide shells. This is solely a
comparison normalization, not a data/model correction. Histogram colour step
.05mag; reported medians are interpolated estimates, not millimag precision.

| z shell | DESI median g−r SGC/NGC | Mock median at matched z SGC/NGC |
|---|---|---|
| .15–.25 | .975 / .989 | .979 / .983 |
| .25–.35 | 1.286 / 1.298 | 1.304 / 1.306 |
| .35–.45 | 1.589 / 1.588 | 1.561 / 1.559 |
| .45–.55 | **1.726 / 1.715** | **1.555 / 1.553** |

The high-shell difference persists after matching fine-redshift composition and
a common numeric r limit: roughly **.17/.16mag redder DESI**. It strengthens the
case for investigating photometric/passband, colour/K-correction and luminosity
modelling. It does not isolate which is wrong, because their passband definitions
are themselves part of the unresolved equivalence question. Remaining within-bin
redshift and magnitude-distribution differences are not removed by this diagnostic.

The published Abacus BGS method uses SDSS/GAMA luminosity-function constraints,
GAMA evolution parameters P=1.8,Q=.7, and colour-dependent GAMA K-corrections;
it is not automatically a fit to this Loa population. This is a concrete
recipe to investigate, **not a verified identification of the exact v0.1 files**.
[Smith et al. 2024, sections5.1–5.2](https://academic.oup.com/mnras/article/532/1/903/7695304).
The newer [DESI DR2 luminosity-function study](https://arxiv.org/abs/2511.01803)
reports DESI-derived K-corrections superseding GAMA prescriptions. Neither paper
alone proves the cause of our discrepancy or authorizes a recipe substitution.

## Required closure work

1. Recover the original targeting-catalogue join/provenance for Loa TARGETIDs,
   including Gaia and REF_CAT/SGA, and pin the actual targeting revision. Replay
   selection including recovery branches; distinguish updated LSS photometry
   from targeting-time photometry.
2. Pin the official mock forFA/altMTL/full-HPmap production commands, map files,
   bad-fibre/tile/exposure masks and the raw CutSky generator recipe. Establish
   R_MAG_APP/G_R_OBS passbands, luminosity/colour evolution and K-corrections.
3. Compare parent and selected populations conditional on fine z, PHOTSYS,
   magnitude, colour and fibre magnitude where physically defined. Reproduce
   stage-by-stage retention. Missing mock observables require an explicit
   observation model; identical column names/cut numbers are not enough.
4. Repair only demonstrated selection/implementation/model defects. Then rebuild
   affected training inputs and validate the impact on encoder/OOF/posterior
   calibration before a versioned replacement. No empirical count adjustment,
   retraining or VAC mutation performed in this audit.

The actual generator and historical selection crosswalk are still unresolved;
full selection alignment and science release remain unqualified. The measurable
checks above completed successfully; no additional phase truth was opened.
