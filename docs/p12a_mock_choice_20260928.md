# Mock choice after photometry/LF diagnostics — 2026-09-28

## Decision

User authorizes release of the existing VAC labelled provisional with explicit
high-redshift inference warnings until improved mocks resolve the issue. This
supersedes withholding the release pending mock replacement; it does not turn
mock calibration into demonstrated real-data coverage. No publication action
or catalogue modification occurred in this diagnostic task.

Keep the frozen P12-A / historical training population for that existing
provisional product, preserving its validation lineage. Canonical Abacus v1 is
the leading replacement-development parent, NOT a validated replacement today.
No tested population presently establishes full .15–.55 Loa statistical and
galaxy–matter parity. Do not alter Loa to reproduce a deficient mock n(z).

## New controlled experiment

`workflows/p12a_vac/parent_evolution_test.py` ran on CPU allocation59018299,
nid004210, completed exit0, step20s. Deterministic stride20, exposed ph000 only,
common NSIDE256 mask and numerical12<=r<19.5. Stored v1 apparent r is changed
only by -(Q-.67)(z-.1). All parent galaxies are eligible, including fainter
rows; no changes persisted to inputs. Stored rest colour split .75 assumes the
published reference convention. No HOD/LF/colour refit or observation simulation.

Selected-count ratio to unchanged v1 (photometric South / North):

| z | Global Q=.78 | Colour Qred=.23, Qblue=1.59 | Q=0 control |
|---|---:|---:|---:|
| .15–.25 | 1.010 / 1.010 | 1.033 / 1.031 | .942 / .938 |
| .25–.35 | 1.032 / 1.036 | 1.071 / 1.063 | .818 / .810 |
| .35–.45 | 1.077 / 1.083 | 1.227 / 1.278 | .619 / .596 |
| .45–.55 | 1.189 / 1.186 | 2.218 / 2.368 | .325 / .305 |

Baseline counts exactly reproduce the independent earlier joint histogram for
both old and v1 parents. Source size/mtime unchanged; script hash, job and count
arrays recorded in evidence/p12a_evolution_test_20260928. No significance claim
from a deterministic sample of one phase. Availability of fainter rows is
recorded, but completeness of all newly eligible faint galaxies is unproven.

Published Q values: https://arxiv.org/html/2511.01803, section3.2. These are
sensitivity interventions, not self-consistent replacements for a population
calibrated with another evolution prescription. Global Q makes v1 surplus
worse; colour-dependent Q alone makes it much worse. Q=0 suppresses the tail
strongly but has no independent justification as a correction.

## Synthesis of tests

- Old raw parents are already deficient versus completeness-corrected Loa at
  high z (roughly .55–.59 in .45–.55). Downstream losses cannot create support.
- Prior fixed-row K/E mapping swap increased old high-z counts ~2.5–2.8 times
  and reduced v1 with the old mapping ~59–62%. This demonstrates sensitivity,
  not attribution to K alone: M and colour conventions differ between parents.
- v1 counts exceed corrected Loa even before assignment, but conditional colours
  are not uniformly improved. High-z residuals persist in common z/r cells.
- Recovered input LF has a sharper extreme bright-end decline than independent
  published DESI LF (see p12a_published_lf_20260928.md). These input tables are
  not a measurement of the realised parent LF; its closure is still open.
- Ashley Ross's reply supplied by the user identifies the input Abacus SecondGen
  high-z problem as known and unresolved. This supports upstream investigation,
  but neither identifies the exact LF/K cause nor validates a DR3 replacement.

## Candidate assessment

| Population | Current use | Reason / remaining limitation |
|---|---|---|
| Historical training mocks | Existing provisional VAC and reproducible baseline | Frozen model calibration exists; high-z transfer bias unresolved. |
| Canonical Abacus v1 | Preferred replacement candidate to validate | Full range support; replayed DESI photometric mapping; surplus and colour mismatch remain; HOD branch and observation/truth closure open. |
| Internal Uchuu Y3-v2.0 tested products | Independent comparison below their support limit | Previously measured zero common-sky rows at .50–.55; not a full-range replacement; host/truth mapping unverified. |
| Tested GLAM/HOLI products | Observation/count benchmark | Tested GLAM tail cutoff; HOLI n(z) calibration and placeholder magnitudes do not establish environment-training parity. |
| Future DR3 mocks | Reassess when an identified product is validated | No tested full-range fix established by the supplied communication. |

## Loa treatment and release requirements

Apply additional Loa cuts only to correct a demonstrated error or define a
scientifically intended sample, with equivalent mock selection and renewed
validation. Thinning/reweighting observations to fit mock n(z) would redefine
the tracer population and generally change inferred environments. Completeness
weights are diagnostic population estimates, not substitute physical neighbours
for the encoder. Retain official catalogue selection and documented known
limitations unless evidence independently requires correction.

For provisional release, document training-mock version, selection, redshift
range, and known high-z mismatch in metadata and README; provide redshift-resolved
quality/support information. A warning at z>=.35 can identify the observed onset,
with particularly severe .45–.55 discrepancy; it must not imply that lower-z
real-data coverage is proven. Keep the original posterior values; no ad hoc
probability repair. Replacement requires independently validated population,
observation model, truth mapping and refreshed posterior calibration.

Remaining causal work: realised parent LF in complete cells; independent Loa
K/E luminosities; separate K and evolution interventions with complete faint
support; selected-stage comparison from the pinned forward pipeline; clustering
and galaxy–matter validation. These remain open, not concealed by the choice of
v1 as the leading candidate. No claim that the photometry/LF investigation is
fully closed.
