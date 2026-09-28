# P12-A VAC implementation handoff (2026-09-24)

Status pointer (2026-09-28): canary and full-survey technical products exist
and remain provisional. Use `p12a_vac_operations.md` for CLIs, FITS columns,
and quality bits. This file is the chronological log and is not rewritten
when later notes supersede a pending step.

Active science plan: `../../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md`.
Alignment subplan: `plan_mock_to_loa_alignment_investigation_20260926.md`.

## Implemented: bounded V0 preflight

`workflows/catalog/p12a_vac_preflight.py` verifies the eight small frozen
candidate artifacts by size and SHA256 without deserializing checkpoints.
It binds the candidate/source/phase-registry records, checks the seven-feature
order against the verified FMPE completion marker, inventories archived Loa
catalogue metadata and produces a conservative phase-use ledger. Candidate
replication phases are not automatically authorized or fresh blind data.

No FITS, particle, target array, reserved confirmation payload or large source
catalogue is opened. Recorded catalogue hashes are explicitly not live verification.
An artifact match cannot pass the coordinate, response, golden-mock or release gates.

The live run is archived as `docs/p12a_vac_preflight_report_20260924.json`:
all eight small frozen artifacts verified; `ready_for_desi_canary` remains false.
Codex's sandbox process launcher was missing its temporary
`codex-linux-sandbox` executable, so this one read-only run used the reviewed
unsandboxed command path. This runner problem is separate from Slurm and from
the scientific gates below.

Run in cosmic_env from GraphWeb_DESI (small metadata/checkpoint hashing only):

```bash
/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python -m unittest discover -s tests -p test_p12a_vac_preflight.py
/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python workflows/catalog/p12a_vac_preflight.py --illustris-root ../TNG/Illustris --output p12a_preflight_20260924.json
```

The output must not exist; reports are never overwritten. Exit 2 is expected
while handoff gates remain pending, including if small artifact checks all pass.
Other exceptions are execution failures, not a scientific result. Archive a
completed small report with the run/source identity after inspecting it.

## V1 schema crosswalk (implementation contract; not field qualification)

| Required item | Source/operation | Constraint |
|---|---|---|
| Galaxy identity | TARGETID from Loa full HPmapcut data | Unique output identity; no row-number identity across releases |
| Angular geometry | RA/DEC | Preserve units and target association |
| Observed redshift | Full-data Z_not4clus candidate, with ZWARN and redshift-quality fields | Freeze real-data success/quality policy before selection; do not copy the mocks' all-success assumption |
| Count input | Successful selected BGS BRIGHT positions on frozen lattice | No automatic PIP/ZFAIL/SYS weighting |
| Angular support | Full HPmapcut randoms and imaging masks | Same response construction as frozen training; holes remain holes |
| Expected counts / ntilde | Frozen training selection and coordinate-volume contract | Clustering-random Z is forbidden as a selection fit source |
| Completeness / redshift success | FRAC_TLOBS_TILES, FRACZ_TILELOCID, mod_success_rate and quality metadata as available | Audit semantics; not new model channels without retraining/validation |
| Cartesian radius | Frozen P12 Planck18 convention pending V0-C outcome | Do not mix native/E2E coordinates with old selection or support distances |
| Encoder output | Frozen epoch-20 U-PATCH, same context and normalization | Resolve exact checkpoint/transform ancestry and repeat final-checkpoint parity |
| Posterior conditioning | Three base predictions, z, log ntilde, cap, log1p boundary distance | Exact seven-feature order; units Mpc and Mpc^-3 |

Source authority: Illustris `configs/p10_response_sources_v1.json`. It names
Loa DA2/LSS/loa-v1/LSScats/v2.1: full data/randoms under the root and clustering
products under PIP. Do not append `/PIP` indiscriminately to full-catalogue paths.
This historical inventory must be re-frozen for P13. A weight column's presence
is not proof that its semantics match training.

## Remaining execution sequence

1. Run/inspect preflight; resolve missing or changed artifacts without modifying
   frozen candidates. Complete full encoder/OOF transform ancestry.
2. Bind licensed training-phase row samples and native host-label sources for
   the coordinate audit; establish criteria and bounded compute request.
   The first ph002 audit is implemented in Illustris; see
   `../TNG/Illustris/docs/p12a_coordinate_sample_run_20260924.md` from this
   repository root. Nine synthetic tests pass; approved job58823054 completed
   sampled coordinate/native-label checks successfully on ph002-005. Results:
   `../TNG/Illustris/docs/p12a_coordinate_sample_results_20260924.md`.
   The imported ph000, lineage/response and Loa checks remain open. All VAC
   compute is now authorized; these samples do not yet qualify the Loa adapter.
3. Implement full data/response adapter only against the resolved coordinate
   and success-selection contract; test IDs, counts, support and interpolation.
4. Replay a licensed golden mock through the exact adapter/checkpoints, then a
   bounded Loa trial. No ph001 reuse for tuning or new blind claims.
5. Additional-phase replication/misspecification before public science release.

No Slurm submission is included. The coordinate audit cannot be replaced by the
presence of matching Planck18 strings. Until its outcome is evidenced, no real
DESI predictions or corrected-model/retraining decision is claimed.

## Live Loa source check completed (2026-09-24)

`workflows/catalog/p12a_loa_source_audit.py` ran in CPU job58823514.
`docs/p12a_loa_source_audit_20260924.json` records matching full/clustering data
SHA256 and matching sizes/row counts/ordered columns for all36 random files.
The random-file content hashes and release success policy remain pending.
This is a source/schema audit, not permission to launch inference.

The installed LSS BGS success helper uses ZWARN==0 and DELTACHI2>40, whereas
our historical `build_bgs_maglim_catalog.py` uses >=25 plus GALAXY. Freeze the
Loa release policy explicitly; do not reuse that legacy selection silently.
The updated Illustris handoff results also close the ph000 import check and
verify embedded checkpoint transforms. Population and response replay gates
remain open; see `../TNG/Illustris/docs/p12a_handoff_followup_results_20260924.md`.

## Actual quality-cut census and replay blocker (2026-09-24)

The full census confirms that the historical catalogue satisfies its own
ZWARN0/DELTACHI2>=25/GALAXY cuts. The installed LSS >40 helper is different,
but that alone is not a defect or a retraining trigger. Quantified threshold
and spectral-type effects: `p12a_quality_cut_results_20260924.md`.

The current scientific blocker is now the Illustris full-fit halo24 encoder's
failed context-growth/subdivision gates. Matching original precision restores
exact stored OOF predictions, so source replay is not the blocker. A frozen-
weight halo48 control passes dense/sparse/edge checks but is not calibrated for
production yet. Follow the active VAC plan's larger-context summary/posterior
revalidation path before declaring the DESI adapter ready. Preserve the Loa
source crosswalk and successful response checks; do not restart E2E work.


2026-09-25: conditional coverage and eight-core observer golden replay pass;
Loa full-data/random content hashes verified. Supersedes pending replay status.
Next is bounded Loa input construction/QA then diagnostic inference; no further
retraining indicated. Full release and independent confirmation remain open.
See Illustris `docs/p12a_golden_conditional_results_20260925.md` and the explicit
bounded-trial handoff manifest. That was the pre-trial state; the completion update below supersedes it.


## 2026-09-25: real Loa canary VAC completed

Supersedes earlier pending coordinate/replay/inference status above. Frozen
halo48 posterior applied to16 preselected Loa cores;5,615 unique TARGETIDs,
5,602 supported posteriors with512 joint draws,13 flagged/null. Independent
serialized-product QA passes. Existing full-data/random catalogues reused;
new P12 field/response inputs built. See `p12a_loa_canary_20260925.md` and
`evidence/p12a_loa_canary_20260925/`. V5 full-footprint scale-out remains next,
after distribution/selection closure and representative shard benchmarking.
No new permission required for authorized production compute. Science release
remains gated by replication/robustness; no confirmation phase was accessed.


2026-09-25 full-survey continuation: input/cut audit found an unresolved high-z
count excess, not an accidental extra cut in our processed mocks. See
`p12a_selection_distribution_audit_20260925.md`. Restartable core checkpoints
and a128-shard ownership plan implemented; benchmark and exact restart pass.
Full-footprint output remains provisional, with explicit selection-shift flag.
