# NERSC / Perlmutter (DESI) agent instructions

You are helping with DESI- and NERSC-centric work on Perlmutter.

## Operating principles

- Prefer **safe, reproducible** guidance: show commands, expected outputs, and failure modes.
- Assume **compute nodes may not have internet**. Avoid plans that require `pip install` during a running job.
- Be explicit about **Slurm details** that commonly break jobs at NERSC: `--account`, QOS selection, CPU vs GPU repos, and using `srun` inside allocations.
- Default to using the **right filesystem** for the job:
  - Use CFS for durable project/user storage.
  - Use PSCRATCH for large temporary I/O; warn about purge policies.
- Avoid expensive metadata operations (e.g. `find` over million-file trees). Prefer targeted globbing, `tar`, or curated path lists.

## When the user asks “run this on Perlmutter”

- First determine context: login vs compute, CPU vs GPU, interactive vs batch, and what account/repo to use.
- Provide an `salloc`/`srun` or `sbatch` path, plus a minimal `job.sbatch` when appropriate.
- Call out the “top gotchas” early (account suffix, QOS namespaces, `srun` requirement, no-internet on compute).

## DESI conventions

- When referencing productions/releases, prefer a **script or authoritative path** rather than hardcoding names (they drift).
- Use canonical storage conventions where possible:
  - User outputs: `/global/cfs/cdirs/desi/users/$USER/`

