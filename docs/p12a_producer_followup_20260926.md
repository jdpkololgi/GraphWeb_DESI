# Producer Abacus and Holi follow-up

Continues the Mock-to-LOA Alignment Investigation, W1/W2. Only already exposed
Abacus ph000 and previously inspected Holi seed0000 were opened. No catalogue
repair, new training, fresh confirmation exposure or VAC replacement.

## What can proceed without user input

Readable producer files provide another independent *pipeline-version* check,
not another independent cosmological phase. Their separate count censuses use
the frozen common occupied NSIDE256 sky, galactic SGC/NGC, 12<=R_MAG_APP<19.5,
and .15<=RSD redshift<.55. Numerical photometry cuts are not proof of physical
passband parity; raw/forFA populations and observed Loa are different stages.
Neither branch may be added to the other without knowing its intended use.

## Producer identity and mapping

Root: `/global/cfs/cdirs/desi/users/jpiat/abacus_mocks/cutsky/z0.200/AbacusSummit_base_c000/ph000`.
N/S raw rows: 71,857,412 / 75,711,865. Y3 forFA rows:
25,023,399 / 26,371,919. Raw N has `v_rel`, S has `vel_disp`;
headers do not pin generation parameters. These products differ from canonical
v1, which has 74,748,522 rows.

The 2,048 evenly spaced row probes per file found all sampled forFA BGS_TARGET
values zero. Full flag census is reported below when completed. A numerical
r<19.5 selection here does not assert that the saved BRIGHT bit is populated.
Old-table colour residual RMS on selected raw samples is 0.054/0.081mag (N/S),
confirming that the v0.1 GAMA mapping does not reproduce these products.

Within the first 10,000 raw rows per branch, 6,175 CEN=1 rows share exact
HALO_ID/RA/DEC/Z keys. Their S-minus-N apparent-r differences have RMS0.343mag,
absolute-r RMS0.403mag, and rest-colour RMS0.259mag. This is a bounded shared
central-location comparison, not a representative population estimate or proof
of galaxy-identity semantics. It establishes overlapping locations with different
assigned properties, ruling out treating N/S as disjoint geographic partitions.
It does not establish which random seeds, passbands or HOD/LF recipes differ.

Evidence: [probe](evidence/p12a_alignment_execution_20260926/PRODUCER_PROBE.json),
[central locations](evidence/p12a_alignment_execution_20260926/PRODUCER_CENTRAL_PAIRS.json).
No matching `jpiat` or `galaxy_cut_sky_N/S` references were found in installed
LSS/main/scripts/mock_tools Python/shell files. The producer root and `misc`
listing contained data products, not the missing source/configuration.

## Readable Holi parents and source-level mechanism

`cai/holi/webjax_v4.82/seed0000` contains readable BGS (6,074,657 rows) and
BGS-NONKP (43,172,697 rows) parents. Both have only RA/DEC/Z/Z_COSMO/NX, without
physical photometry or halo IDs. This supersedes any broad statement that no
readable Holi parent exists; the delivered forFA0 access restriction is separate.

The archived installed `prepare_holi_bgs.py` explicitly names this DA2 root and
these filenames. It thins NONKP to Loa ANY n(z), then to Loa BRIGHT n(z). It
separately augments the BGS parent with an ANY-tail sample at .5<z<.61 and
constructs the BRIGHT-21.35 cosmology subsample. Magnitudes -22/-21/-19 are
assigned as category constants. Close n(z) can therefore be a construction
criterion, rather than independent evidence that the physical luminosity model
matches Loa. This is a source-level mechanism, not a claim about the executed v2.

A specific provenance mismatch remains: the inspected script defaults to DA3
and loops over seeds3–99, excluding seed0000; the inspected pipeline wrapper
names `holi_bgs_bmask`, not the delivered `holi_bgs_v2` product. We cannot treat
these sources as an execution record for our mock0. No downloaded or shared
production script was executed; source was inspected read-only.

Evidence: [source audit](evidence/p12a_alignment_execution_20260926/HOLI_SOURCE_FOLLOWUP.json).

## Remaining inputs to physical regeneration

The producer's executed source/configuration, LF and HOD tables, K/E/passband
conventions, random seeds, and N/S combination rule are needed before controlled
P/Q or photometric experiments can be interpreted as changes to this generator.
The Holi delivered-product link likewise needs its actual preparation version.
These are provenance gaps, not requests for more compute approval. Common-sky
screening can continue without them; production qualification cannot infer them
from matching counts alone.

## Canonical v1 crosswalk

An initial exact HALO_ID/RA/DEC/Z prefix crosswalk found zero matches. Relaxing
only the RSD-Z key resolves 6,322/6,260 shared central angular locations for
producer N/S versus canonicalv1. Their Z_COSMO values are exactly equal, but
RSD-Z differences have RMS0.001284/0.001311 and velocity-component differences
have RMS263–270 in the stored units. Apparent-r differences have RMS0.339/0.334mag.
This is a generation/velocity-assignment difference, not just applying different
flags to an identical galaxy table. The prefix sample cannot measure whole-survey
impact or determine the intended velocity prescription. Keep both conventions
explicit and require the producer's recipe. Evidence:
[canonical crosswalk](evidence/p12a_alignment_execution_20260926/PRODUCER_V1_PAIRS.json).

## Completed census and validation

CPU58908840/nid004153 completed all four files. Counts below are numerical
r-selected common-sky counts; caps in each cell are SGC / NGC.

| Product | .15–.25 | .25–.35 | .35–.45 | .45–.55 |
|---|---:|---:|---:|---:|
| Raw N | 1,227,514 / 2,768,855 | 750,017 / 1,694,035 | 294,615 / 680,812 | 72,558 / 165,865 |
| Raw S | 1,288,700 / 2,905,829 | 817,587 / 1,843,899 | 326,570 / 756,834 | 79,596 / 181,884 |
| Y3 N | 1,218,050 / 2,756,251 | 744,071 / 1,686,343 | 292,265 / 677,880 | 71,974 / 165,147 |
| Y3 S | 1,278,766 / 2,892,628 | 810,978 / 1,835,656 | 324,014 / 753,518 | 78,940 / 181,110 |

All25,023,399 N-Y3 and26,371,919 S-Y3 rows have BGS_TARGET=0. This
does not establish an error in the historical P12 inputs: these are distinct
producer products. Direct BRIGHT-bit selection on these files would select
nothing; numerical-r counts are explicitly a diagnostic, not a targeting replay.

Raw high-shell ratios to original v0.1 are2.490/2.524 (N branch) and
2.732/2.768 (S branch), for galactic SGC/NGC. These parents have additional
bright-tail support; whether their observed populations match Loa is not tested.

[Comparison figure](figures/p12a_alignment_execution_20260926/producer_branches.png).
All census histogram sums, shared-mask hashes, source/helper hashes and probe
hashes passed. Four relevant existing unit tests passed. Plot visually inspected;
matched row-wise axis scales used. No statistical qualification or physical
recipe identification is implied. See `PRODUCER_VALIDATION.json` and plot receipt.

### Public producer context and focused clarification (2026-09-26)

User identifies directory owner as Jade Piat. "Producer mocks" was our local
shorthand for files in her users/jpiat/abacus_mocks tree, not an official release
name. Her public COLOURS talk (10 June2025), PDF pages18–22 and33–34, describes
AbacusSummit, a magnitude-dependent HOD, a z=0.2 snapshot with magnitude evolution,
and velocity-evolution work. This supports a specific question about velocity
prescriptions; it does not identify the inspected June2026 files or explain their
Loa agreement. Source:
https://indico.ijclab.in2p3.fr/event/11110/contributions/37832/attachments/25782/37994/PIAT_Relativistic_effects_DESI_COLOURS.pdf

Do not conflate bright/faint tracer splits in the dipole study with DESI targeting
classes BGS_BRIGHT/BGS_FAINT. No public exact-file generation/Loa selection recipe
located in this targeted search. Web PDF text available; screenshot fetch failed.
No claims inferred from unread figure contents.

Clarification priority: recommended existing full-BRIGHT Loa product; exact
N/S/passband and magnitude/footprint/target/quality selections; any input n(z),
LF/evolution/K-correction changes; observation processing and validation domain;
velocity convention and halo mapping. A short documented procedure, parameter
values and catalogue/version pointers can unlock selection replay: obtaining an
entire code repository is not a prerequisite merely to test existing products.
Our parent-count screen establishes extra bright-tail support, not statistical
Loa agreement. Selecting a few flags is a hypothesis, not a demonstrated solution.
Slack message drafted for the user; no external message sent.
