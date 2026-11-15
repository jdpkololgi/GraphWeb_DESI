# GraphWeb_DESI Runbook

This runbook lists validated workflow entrypoints and launch commands.

## Environment setup

Recommended shell setup for GraphWeb workflows:

```bash
source ~/.bashrc
conda activate cosmic_env
desienv
```

If `desimodel` import remains unavailable in your shell context, use the
canonical script paths and environment used in validation notes. The catalog
workflow also includes a hardcoded `desimodel` fallback path.

## Canonical workflows

### Catalog assembly

Canonical:

```bash
python workflows/catalog/load_catalog.py --help
```

Compatibility wrapper:

```bash
python load_catalog.py --help
```

Example run:

```bash
python workflows/catalog/load_catalog.py --out-dir .
```

### Graph inference pipeline

Canonical:

```bash
python workflows/graph_inference/graph_catalog.py --help
```

Compatibility wrapper:

```bash
python graph_catalog.py --help
```

Example no-plot cache-only probe:

```bash
python workflows/graph_inference/graph_catalog.py --cache-mode cache-only --no-summary-plot
```

If model checkpoint path differs from config default:

```bash
python workflows/graph_inference/graph_catalog.py --model-path /path/to/trained_gat_model.pth
```

### Gudhi/cuGraph/Jraph wedge inference

This path is for DESI wedge eigenvalue regression with an Abacus-trained Jraph
model. It is not a replacement command for the all-sky GAT VAC workflow above.
Run it on Perlmutter from an environment with DESI, Gudhi/cuGraph where needed,
JAX/Haiku/Jraph, and access to the Illustris training repo.

Build the bright BGS catalog if it does not already exist:

```bash
python workflows/catalog/build_bgs_maglim_catalog.py \
  --out-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits
```

Build the full Mpc-parity graph and export graph features:

```bash
python workflows/graph_construction/build_desi_bgs_gudhi_graph.py \
  --catalog-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc \
  --alpha-sq inf \
  --coord-units mpc

python workflows/graph_construction/desi_graph_features_cugraph.py \
  --metadata-path /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_delaunay_metadata.json \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc \
  --out-prefix desi_bgs_cugraph
```

Subset the expanded Abacus-parity DESI wedge and run inference:

```bash
source workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt

python workflows/graph_construction/subset_desi_graph_wedge.py \
  --graph-metadata /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_delaunay_metadata.json \
  --catalog-path /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --parent-gnn-arrays /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_bgs_cugraph_gnn_arrays.npz \
  --parent-gnn-metadata /pscratch/sd/d/dkololgi/graphweb_desi/outputs/gudhi_hemi_full_alphasq_inf_seed42_bright_mpc/desi_bgs_cugraph_gnn_metadata.json \
  --out-dir /pscratch/sd/d/dkololgi/graphweb_desi/outputs/desi_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc_from_fullgraph \
  --out-prefix desi_delaunay_wedge_expanded_ra120_160_dec14p5_30p6_z0p2_0p3_bright_mpc \
  --ra-min 120 --ra-max 160 --dec-min 14.5 --dec-max 30.6 --z-min 0.2 --z-max 0.3

python workflows/jraph_inference/jraph_infer_desi_wedge_from_gnn_npz.py \
  --abacus-run-dir "${ABACUS_RUN_DIR}" \
  --calibration-cache "${CALIBRATION_CACHE}" \
  --abacus-gnn-arrays "${ABACUS_WEDGE_GNN_ARRAYS}" \
  --desi-gnn-arrays "${DESI_GNN_ARRAYS}" \
  --desi-gnn-metadata "${DESI_GNN_METADATA}" \
  --desi-global-node-ids "${DESI_GLOBAL_NODE_IDS}" \
  --desi-wedge-catalog-npz "${DESI_WEDGE_CATALOG_NPZ}" \
  --output-dir "$(dirname "${INFER_DIR}")" \
  --run-name "${INFER_RUN_NAME}"
```

`INFER_DIR` in `workflows/catalog/JRAPH_INPUTS_expanded_wedge.txt` is the full
intended run directory, so the command passes its parent as `--output-dir` and
the basename as `--run-name`. Jraph inference uses `ILLUSTRIS_ROOT` to find
`shared/graph_net_models.py`; the GAT classifier uses `ILLUSTRIS_REPO_ROOT` from
`shared/config_paths.py`.

### DESI absolute-magnitude support diagnostic (experimental)

`workflows/sbi_inference/desi_absmag_kcorr.py` evaluates whether a luminosity
feature is worth further parity work. It joins BGS photometry by `TARGETID`,
calls the official DESI LSS `add_dered_flux` and `add_ke` implementations, and
writes k+e-corrected `ABSMAG_RP1`. It then prints marginal magnitude and colour
statistics beside a fixed raw Abacus cut-sky sample.

Run this in the DESI environment, not `cosmic_env`:

```bash
source /global/common/software/desi/desi_environment.sh main

python workflows/sbi_inference/desi_absmag_kcorr.py \
  --wedge /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_maglim_bright_galaxy_zwarn0_dchi2ge25.fits \
  --full /global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/LSScats/v2.1/BGS_BRIGHT_full.dat.fits \
  --out /pscratch/sd/d/dkololgi/graphweb_desi/catalogs/bgs_absmag_rp1_gate_sub.fits \
  --zlo 0.15 --zhi 0.55 \
  --max-rows 400000 --seed 42
```

The `--lsscode` checkout must contain both `LSS/py/LSS` and `DESI_ke`; its
default is the complete shared NERSC checkout used when this workflow was
developed. The input named `--wedge` must provide `TARGETID`, `Z`,
`TARGET_RA`, and `TARGET_DEC`; `--full` supplies `FLUX_G`, `FLUX_R`, and the
corresponding Milky Way transmissions. The current script does not apply an
RA/Dec cut to the DESI input despite the historical option name.

Operational constraints:

- The default 400,000-row sample is selected reproducibly after the positive
  flux and half-open `zlo <= Z < zhi` cuts. `--max-rows 0` processes every
  selected row; the 8.7-million-row parent was projected to take about four
  hours in the single-threaded `add_ke` solve.
- The output is a support-study subsample, not a deployment catalogue.
- The FITS output is written before the Abacus comparison runs. A later
  comparison failure does not invalidate the k+e columns, but it does mean the
  diagnostic did not complete.
- The Abacus path and its fixed sky/magnitude selection are currently
  hard-coded. There is no CLI pass threshold or nonzero exit status for a
  support mismatch.

Interpret the printed offsets as exploratory diagnostics only. The present
comparison uses different DESI and Abacus footprints and unmatched redshift
distributions, and it does not prove that DESI `ABSMAG_RP1` and Abacus
`R_MAG_ABS` share the same band, reference-redshift, distance-modulus, and
evolution-correction conventions. Before adding luminosity to production
inference, compare the exact training and inference selections in redshift
bins (or with matched `n(z)`), establish the magnitude convention explicitly,
and define quantitative acceptance thresholds.

### Utility workflows

Canonical:

```bash
python workflows/utilities/galaxy_catalog.py --help
python workflows/utilities/investigate_edges.py --help
```

Compatibility wrappers:

```bash
python galaxy_catalog.py --help
python investigate_edges.py --help
```

## Notes

- Active workflows and statuses: `ACTIVE_WORKFLOWS.md`
- Reorg notes: `docs/migration/WORKFLOW_REORG.md`
- Shim policy: `scripts/SHIM_DEPRECATION.md`
- Prefer canonical `workflows/...` paths in new docs and launch scripts.
