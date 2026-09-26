# Reproduction and execution

Working directory: GraphWeb_DESI. Python:
`/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python`.
Clear PYTHONPATH/PYTHONHOME/PYTHONUSERBASE/LD_PRELOAD; set PYTHONNOUSERSITE=1.
Source and helper hashes are in completion receipts; CFS read-only.

CPU interactive58921665/nid004148,1h request,64 allocated CPUs. DR1 and boundary
steps used8 CPUs each, overlapped within the same node with sufficient resources.
The allocation shell retained the foreground DR1 step, avoiding idle-shell expiry.
Unrelated E2E GPU58920768 left untouched. Released after both completed and plots
rendered. Scheduler/steps all COMPLETED0; see COMPUTE.json.

Commands (prepend srun -N1 -n1 -c8 --cpu-bind=cores inside the allocation):

```
python -u workflows/p12a_vac/alignment_dr1_controls.py --root /pscratch/sd/d/dkololgi/graphweb_desi/outputs/alignment_dr1_controls_20260926
python -u workflows/p12a_vac/alignment_boundary_controls.py --root /pscratch/sd/d/dkololgi/graphweb_desi/outputs/alignment_boundary_controls_20260926
```

Second used srun --jobid=58921665 --overlap from the login shell. Separate input
and output paths; both immutable source scripts throughout execution. Completed
roots fail closed on rerun. Fresh output root required for reproduction.
Compact JSON/NPZ outputs copied into dr1/ and boundary/ here.

Bounded10k-row prefix diagnostics (lightweight login checks):

```
python workflows/p12a_vac/alignment_velocity_control.py --root docs/evidence/p12a_independent_controls_20260926
python workflows/p12a_vac/alignment_rsd_speed_check.py
python -m unittest discover -s tests -p test_alignment_controls.py
python workflows/p12a_vac/plot_alignment_controls.py --dr1 docs/evidence/p12a_independent_controls_20260926/dr1 --boundary docs/evidence/p12a_independent_controls_20260926/boundary --out docs/figures/p12a_independent_controls_20260926
```

Five focused tests pass; conservation and provenance validation recorded in
VALIDATION.json. Final compact-data plot rerender on login only clarified that
the original Loa/mock mask includes non-DR1 sky. Both PNGs visually inspected.
Graphify first attempted from home (no graph), then correctly run in GraphWeb
using cosmic_env. Final AST update/global refresh performed. No scientific run
failed, no retries, new phase exposures, catalogue repairs or model fitting.
