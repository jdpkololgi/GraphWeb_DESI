# Research context — GraphWeb_DESI

## Role

GraphWeb_DESI is the DESI observational application layer of GraphWeb. It
processes DESI BGS catalogues, constructs survey-aware graphs, applies
simulation-trained models, and produces inferred eigenvalue and cosmic-web
data products.

The simulation and methods source of truth is the sibling `Illustris`
repository. Its `SCIENCE_LOG.md` is the shared live scientific record.

## Scientific goal

Apply calibrated inference of the ordered tidal-tensor eigenvalues

`lambda_1 <= lambda_2 <= lambda_3`

to DESI BGS galaxies, then use the resulting Value-Added Catalogue for
environmental science. Discrete void/wall/filament/cluster labels are derived
products, not the fundamental inference target.

## Current interface to Illustris

- The current production anchor is the Illustris G3 union-graph NPE model:
  attentional Jraph GraphNetwork plus FlowJAX posterior.
- Training parity uses AbacusSummit HOD cutsky mocks with CACTUS T-web labels,
  7 Mpc/h smoothing, and `lambda_th = 0.2`.
- Preserve the canonical ordered-increment target representation and the
  Abacus-compatible graph feature and coordinate conventions.
- Use `ILLUSTRIS_ROOT` / `ILLUSTRIS_REPO_ROOT` configuration rather than
  embedding machine-specific repository paths.

## Current DESI state

- DESI BGS DR2 Loa / FastSpecFit integration is the observational basis.
- Approximately 112,755 galaxies have zero-shot eigenvalue inference with a
  roughly 99.7% TARGETID match to FastSpecFit.
- Environmental closure tests use observables such as quenched fraction and
  sSFR at controlled stellar mass; CIGALE-HZ rejoin is an optional property
  swap for those plots when FastSpecFit SFRs are unsuitable.
- Abacus-domain transfer gates (GNN vs GBM capacity, RSD/luminosity aperture
  features, mock–DESI n(z)) live under `workflows/sbi_inference/` and inform
  feature work — they are not DESI VAC acceptance criteria.
- TARP/SBC validate the model inside the Abacus domain; DESI closure tests are
  the necessary truth-free observational check.

Operational entrypoints in this repo (see `ACTIVE_WORKFLOWS.md` / `RUNBOOK.md`):

- GAT VAC-style classification: `workflows/graph_inference/graph_catalog.py`
- Jraph wedge point regression:
  `workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py`
  (expanded-wedge wrapper: `workflows/catalog/run_infer_expanded_wedge_mpc.sh`)
- FlowJAX/SBI wedge posteriors: `workflows/sbi_inference/`
  (README there is the detailed parity/runbook companion)
- Gudhi/cuGraph Slurm stages: `workflows/catalog/sbatch_desi_bgs_bright_pipeline.sh`
  (retarget before treating outputs as the expanded Mpc-parity wedge)

Low-z GAT catalogs live on pscratch via `GRAPHWEB_CATALOG_*` in
`shared/config_paths.py`; do not reintroduce large FITS products under home.

## Active scientific direction

The Illustris-side priority is field-level, physics-grounded inference:

`graph encoder -> density grid -> fixed FFT tidal operator -> eigensolver`

This repository should prepare for and evaluate the resulting production
encoder, but should not treat field-level point-estimate results as a DESI
replacement until the associated posterior/calibration gate is passed.

Key caveats:

- The raw field-level tests used a denser-than-DESI wedge. DESI-like number
  density and survey-selection controls are required before transfer claims.
- NPE contains the Abacus training prior; do not stack DESI posteriors for
  population-level eigenvalue PDFs without reweighting or hierarchical SBI.
- Feature-space domain shift and prior mismatch are separate risks.
- T-web is the potential Hessian, not V-web velocity shear.

## Data-product principles

- Preserve TARGETID joins and provenance for every derived prediction.
- Keep coordinate units explicitly documented; the Gudhi graph builder uses
  comoving Mpc by default for Abacus-training parity.
- GAT alpha/Delaunay graphs (Illustris `network` + `graph_catalog.py`) are a
  different construction from the Gudhi/cuGraph Mpc wedge graph. Do not feed
  GAT cache pickles into Jraph/SBI inference.
- Keep observational catalog processing separate from simulation-side model
  training.
- Never overwrite a release/VAC artefact without a reproducible replacement
  and clear provenance.

## Shared workflow

Read the Illustris `SCIENCE_LOG.md` before substantive work. It wins if it
conflicts with this file.

Use the same `[science]` and `[code]` log conventions across both repositories.
Synchronise changes through git, pulling with `git pull --no-rebase` before
pushing and preserving all valid log entries in conflicts.