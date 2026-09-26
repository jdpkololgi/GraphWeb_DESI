# Reproduction and run ledger

Project Python used cosmic_env on allocated CPU nodes. Before Python:

```bash
unset PYTHONPATH PYTHONHOME PYTHONUSERBASE LD_PRELOAD
export PYTHONNOUSERSITE=1
```

Executable: `/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python`.
Working directory: `/global/u2/d/dkololgi/GraphWeb_DESI`.
Allocations used one node,64CPUs,desi account,interactiveQOS,scratch license
and1hour request. COMPUTE_CLOSEOUT.json preserves actual elapsed times,
step exits and nodes. First shell auto-logged out while idle after completed
runs; second allocation explicitly released. Unrelated jobs untouched.

Entrypoints under workflows/p12a_vac take --root (fresh output directory):

| Script | Output root suffix under graphweb_desi/outputs | Additional arguments |
|---|---|---|
| alignment_parent_screen.py | alignment_parent_20260926 | --version v0.1, then --version v1 |
| alignment_version_fingerprint.py | alignment_parent_20260926 | none |
| alignment_uchuu_screen.py | alignment_uchuu_20260926 | none |
| alignment_internal_screen.py | alignment_internal_20260926 | none |
| alignment_release_crosswalk.py | alignment_release_20260926 | none |
| alignment_uchuu_dr2_screen.py | alignment_uchuu_dr2_20260926 | none |
| alignment_uchuu_loa_parity.py | alignment_uchuu_loa_parity_20260926 | --parent internal Uchuu output root |

Scratch prefix: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/`.
Scripts pin source paths and selection contracts. Completed census markers
fail closed. Parent version screen supports verified completed-version skip.
Fingerprint/plot scripts regenerate diagnostic summaries; use separate roots
to preserve a prior experiment.

Plot entrypoints: plot_alignment_screens.py takes --parent (Abacus), --uchuu
(SV3), --internal (GLAM/Holi) and --out. plot_uchuu_dr2_alignment.py takes
--parent (internal Uchuu), --parity (current Loa crosswalk) and --out.
Figures and compact histograms are committed; no raw catalogues copied to Git.

Receipts preserve script SHA256. Initial implementation/results commit9b3316b;
final support checks/plots093cf5a. Match source hashes, not an assumed later
working-tree revision. Upstream snapshots have independent hash manifests.

Failures: Holi target read denied; compute-node Graphify flock524 recovered
on login in cosmic_env; ad-hoc string BGS_TYPE min/max probe corrected.
Host Python3.6 lacks capture_output; ledger writer used compatible subprocess
arguments. Scientific completion uses COMPLETE.json plus reconciled histograms,
not scheduler status alone. VALIDATION.json records checks and limitations.
