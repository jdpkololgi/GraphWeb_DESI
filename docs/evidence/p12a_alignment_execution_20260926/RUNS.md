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


## Producer follow-up: CPU58908840, nid004153

One CPU interactive allocation (1h limit,64 CPUs allocated; census srun used8),
retained PTY with foreground srun to avoid idle-shell expiry. Unrelated GPU
allocation left untouched. Scientific Python cosmic_env, cleared Python/loader
environment; PYTHONNOUSERSITE=1. CFS inputs read-only.

- `python workflows/p12a_vac/alignment_producer_probe.py` (bounded2048 rows/file).
- `python workflows/p12a_vac/alignment_producer_central_pairs.py` (first10k N/S).
- `python workflows/p12a_vac/alignment_producer_v1_pairs.py` (first10k canonicalv1/N/S;
  exact RSD-key dead end, angular-key diagnostic additionally checks Z_COSMO/velocities).
- `srun -N1 -n1 -c8 --cpu-bind=cores /pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python -u workflows/p12a_vac/alignment_producer_census.py --root /pscratch/sd/d/dkololgi/graphweb_desi/outputs/alignment_producer_20260926`
- `plot_producer_alignment.py --root docs/evidence/p12a_alignment_execution_20260926/producer`
  final compact-data rerender uses equal y-scales per row. PNG visually inspected.
- Four relevant unit tests passed; source/helper/mask/histogram hashes verified.

Allocation checker initially used system python3 (3.6), failing before execution
on future annotations; rerunning with cosmic_env Python succeeded. Graphify update
and global refresh used cosmic_env on login. No scientific job failed. Four census
receipts plus histogram evidence archived under producer/. Explicit shell exit
released allocation after census and first plot. No jobs left running by this task.
