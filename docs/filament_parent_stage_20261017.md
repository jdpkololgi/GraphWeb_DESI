# Frozen P12-A filament morphology investigation: first compute gate

Current verdict: **unresolved**. Cached projection statistics do not warrant
retraining, visual relabelling, or a change to the production catalogue.

Read authoritative current SCIENCE_LOG, both repository instructions/contexts,
original A06 manifest, filament audit and slab scripts. Original734 line-like
targets include482 with Pfil<.2 (297wall/185void argmax); these are dependent
projected-tracer statistics, not tidal truth/error rates. Existing slab has no
complete10.34585Mpc sphere. Paired membership/eligibility changes prevent a
slab-thickness conclusion. Do not repeat the prior unpaired interpretation.

## Implementation audit so far

Production source uses strict `draws>.2`, sums three eigenvalue crossings per
joint draw, then averages four memberships into P_WEB. It checks ordered draws
and transforms scaled increments back to physical eigenvalues. Posterior-mean
class is a different quantity. Existing plotting audit indexes P_WEB[:,2], and
cached row/coordinate checks report exact equality. All original rows have
QUALITY32 only. These are code inspection and inherited records, not fresh
full-catalogue/draw validation. Exact tidal operator normalization and host-label
versus galaxy-position semantics still require binding from actual training
artifacts; do not use generic density-Hessian wording as that definition.

Manifest pins the September25 full-survey VAC SHA256 cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc
and the posterior/encoder hashes. Original geometry uses Planck18 redshift-space
Mpc and h=.6766 to convert7Mpc/h. No inference weights or selection changed.

## Prepared extraction (not run)

`workflows/p12a_vac/filament_geometry/prepare_parent.py`:
- Fail closed on full catalogue hash, original saved-row hash, or any A06 model
  binding mismatch; hash catalogue again after extraction.
- Chunk full catalogue, without original RA/z/slab or probability selection;
  retain all rows within320Mpc of (0,1025,0), including unsupported tracers.
- Keep TARGETID, catalogue row, coordinates, original quality and probabilities;
  verify all4459 original targets and positions against saved rows.
- Preserve original inputs; require absent output directory, write READY last.

Registered subsequent geometry: full depths5/10/20/40/80Mpc; rotations0/30/60/90
about the radial-plane y axis; normal offsets-5/0/+5Mpc; half-widths180/225Mpc.
Radii0.5/1/1.5 times7/.6766Mpc; minimum neighbours6/8/12; original shape
threshold.5847730874520835. Parent320Mpc sphere encloses every view plus maximum
radius: sqrt(180^2+225^2+45^2)+15.5188<320. All orientations/offsets fixed before
seeing new neighbours/probabilities. They are sensitivity choices, not exhaustive
views or a calibrated3D classification. Tests pass for coordinate axes, rotation
orthogonality, padding and original slab selection (3tests, no catalogue scan).

This stage only extracts/validates parent; no new morphology results exist yet.
Check whether the VAC-selected parent matches inference context selection before
calling the result a full-context test. Guard radial/sky masks separately from
parent padding; BOUNDARY_MPC alone is not a newly validated3D mask certificate.
Then implement/save fixed-target membership transitions, common cohorts, paired
count-thinning and spatial-block uncertainty. Track original482 suspicious rows
and representative probability-independent controls; never optimize viewing angles.
Fresh core-draw P_WEB reconstruction and training/preprocessing parity remain open.

Mock gate: inventory exposed ph006 controls versus genuinely untouched phases
before reading new payloads. ph006 has been used for selection; it is development,
not blind confirmation. Preserve sealed phases. Need matched matter-derived
truth normalization, host/mass/selection correspondence, mask and RSD controls
before any failure attribution. No retraining experiment is justified yet.

## Requested compute scope

One CPU node,64logical CPUs, accountdesi, QOSinteractive,30minutes maximum,
scratch license, cosmic_env. Parent extraction/source hash verification only;
no GPU, training, new phase exposure or Slurm changes to concurrent jobs.
Release only this allocation when done; retain outputs, no file deletion.
Output proposal:
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_parent_20261017_v1`.

Launch after explicit user authorization, using ordinary host execution only:
`python workflows/p12a_vac/filament_geometry/prepare_parent.py --manifest /global/homes/d/dkololgi/fmsu_tokyo_20261010/a06_desi/manifest.json --out /pscratch/sd/d/dkololgi/graphweb_desi/outputs/filament_parent_20261017_v1`

Active chat still advertises workspace-write/auto-review and restricted profile,
inconsistent with newer user preference for native execution. Do not alter the
profile or bypass security hooks. Graphify query returned a scoped graph plus a
MUNGE socket error; this is not a successful Slurm preflight. No Slurm action taken.

## Execution completed — 2026-10-17

Approved CPU allocation59904097 on nid200020 completed0:0; step9s and total16s.
Ordinary sandbox launch failed before Slurm (bwrap directory permission error);
explicit native execution succeeded without changing profile/security settings.
Allocation released automatically; concurrent jobs untouched.

Parent contains180977 rows and all4459 original targets. Fresh full VAC SHA256
before/after extraction matches cf472cb4e8fb327620629f347115ad26c55a3f985320b293e6753e07f50ebfbc.
Posterior, encoder and model-manifest hashes match A06 bindings. Exact original
ID/RA/DEC/Z/quality/support/Pfil checks passed, geometry within1e-8Mpc. READY
receipt copied to docs/evidence/filament_parent_20261017/READY.json; arrays remain
in isolated Scratch output. No off-slab shape analysis has yet run. Verdict
remains unresolved; parent verification is not evidence that structures are walls.
