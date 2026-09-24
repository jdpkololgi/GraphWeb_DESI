# P12-A VAC implementation handoff (2026-09-24)

Active science plan: `../../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md`.

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
3. Implement full data/response adapter only against the resolved coordinate
   and success-selection contract; test IDs, counts, support and interpolation.
4. Replay a licensed golden mock through the exact adapter/checkpoints, then a
   bounded Loa trial. No ph001 reuse for tuning or new blind claims.
5. Additional-phase replication/misspecification before public science release.

No Slurm submission is included. The coordinate audit cannot be replaced by the
presence of matching Planck18 strings. Until its outcome is evidenced, no real
DESI predictions or corrected-model/retraining decision is claimed.
