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

- The current application baseline is frozen Illustris U-PATCH + P12-A FMPE.
  A 16-core Loa canary and a full-survey provisional VAC exist
  (`docs/p12a_vac_operations.md`). Both set `science_release_ready` false.
  The open scientific gate is mock–Loa population alignment, not a missing
  inference entrypoint. G3/FlowJAX results and the legacy GAT/Jraph entrypoints
  are historical products, not this VAC.
- Active cross-repository plan:
  `../TNG/Illustris/docs/plan_desi_p12a_vac_20260924.md`.
  Population subplan: `docs/plan_mock_to_loa_alignment_investigation_20260926.md`.
  Additional mocks support replication and selection tests, not automatic
  retraining or a replacement catalogue.
- Training parity uses AbacusSummit HOD cutsky mocks with CACTUS T-web labels,
  7 Mpc/h smoothing, and `lambda_th = 0.2`.
- Preserve the canonical ordered-increment target representation and the
  Abacus-compatible graph feature and coordinate conventions.
- Use `ILLUSTRIS_ROOT` / `ILLUSTRIS_REPO_ROOT` configuration rather than
  embedding machine-specific repository paths.

## Current DESI state

- DESI BGS Loa full HPmapcut is the P12-A observational basis
  (`Z_not4clus`, `ZWARN==0`, `DELTACHI2>=25`, `GALAXY`; output
  `0.15<=z<0.55`). That cut is not the GAT low-z FastSpecFit sample and not
  the maglim zall used for Gudhi wedges.
- The provisional full-survey note records 5,436,413 TARGETIDs. In the writer,
  quality bit 32 is set on every row; bit 64 flags `Z>=0.35` for an unresolved
  count excess versus the training mocks. Do not publish it as a qualified VAC.
- The older wedge zero-shot (about 112,755 galaxies, FastSpecFit match near
  99.7%) is the Jraph/FlowJAX sky cut, not the P12-A FITS.
- Environmental closure tests use observables such as quenched fraction and
  sSFR at controlled stellar mass; CIGALE-HZ rejoin is an optional property
  swap for those plots when FastSpecFit SFRs are unsuitable.
- Abacus-domain transfer gates (GNN vs GBM capacity, RSD/luminosity aperture
  features, mock–DESI n(z)) live under `workflows/sbi_inference/` and inform
  feature work — they are not DESI VAC acceptance criteria.
- TARP/SBC validate the model inside the Abacus domain; DESI closure tests are
  the necessary truth-free observational check.

Operational entrypoints in this repo (see `ACTIVE_WORKFLOWS.md` / `RUNBOOK.md`):

- P12-A Loa posterior VAC: `workflows/p12a_vac/` (see `docs/p12a_vac_operations.md`)
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
The bright maglim catalog for Gudhi/Jraph/SBI is a separate zall product
(`TARGET_RA`/`TARGET_DEC`, no redshift or mass cut). Its builder defaults to
the public DR2 zall, not `DESI_ZCAT_FILE`. GAT graphs are 10-d Illustris
network features; wedge Jraph/SBI graphs are 7-d Abacus-style features.

## Active scientific direction

Separate Illustris-side research explores field-level, physics-grounded inference:

`graph encoder -> density grid -> fixed FFT tidal operator -> eigensolver`

The immediate application priority is the calibrated per-galaxy P12-A baseline.
Coherent-field work is not a VAC prerequisite. Whole-survey encoding would not
by itself make independently sampled per-galaxy posteriors spatially joint.
Audit final-checkpoint tiling parity and physical context adequacy separately.

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
  comoving Mpc by default for Abacus-training parity. FlowJAX inference
  aborts on Mpc/h metadata; Jraph inference only warns.
- Keep FastSpecFit `ABSMAG01_SDSS_*` (property join) distinct from LSS
  `ABSMAG_RP1` (`desi_absmag_kcorr.py`); they are different magnitude systems.
- GAT alpha/Delaunay graphs (Illustris `network` + `graph_catalog.py`) are a
  different construction from the Gudhi/cuGraph Mpc wedge graph. Do not feed
  GAT cache pickles into Jraph/SBI inference.
- Jraph visualization (`visualize_desi_wedge_cweb_3d.ipynb`) reads
  `desi_wedge_index_and_preds.npz`. FlowJAX talk visuals (skewer / class-3D
  GIF) read `desi_wedge_flowjax_preds.npz`. Do not swap those products.
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
