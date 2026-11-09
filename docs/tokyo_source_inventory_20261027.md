# Tokyo observational presentation source inventory

These files preserve the October 17–19 presentation work and its original,
frozen numerical inputs. They are historical artifacts: the subsequent
October 19 label-provenance warning and current retraining decisions in the
shared SCIENCE_LOG take precedence. Committing these scripts does not recertify
the old model or its targets. Do not silently rerun them with corrected-label
products substituted for their recorded inputs.

Generated slide exports, arrays, receipts and ZIPs remain outside these
repositories in `/global/homes/d/dkololgi/fmsu_tokyo_20261010/`. They are not
included in these commits. Existing ignore rules exclude generated scientific
arrays/figures in the repositories. No source or scientific record is ignored.

The final corrected package is `tokyo_notation_corrections_20261019.zip` under
that output root. It contains the bold lambda_g notation and the expanded
inferred-tidal-trace axis definition; earlier packages remain preserved.

## Source groups

All scripts are under `workflows/presentation/`:

- `desi_filament_visual_audit_20261017.py` and
  `desi_slab_sensitivity_20261017.py`: original bounded projected-slab audit.
  `docs/filament_visual_audit_20261017.md` and
  `docs/filament_investigation_prompt_20261017.md` preserve that initial
  investigation stage; later morphology reports supersede its open questions.
- `tokyo_wall_and_ssfr_overlay_20261017.py` and
  `tokyo_wall_and_ssfr_overlay_extended_20261017.py`: mock main/backup
  examples and matched continuous/class sSFR summaries. The latter adds the
  measured high-trace bin and uses 59,774 galaxies in 69 common cells;
  the earlier version uses 63,347 galaxies in 74 cells. Preserve this difference.
  Initial report: `docs/tokyo_wall_overlay_20261017.md`.
- `tokyo_galaxy_properties_labels_20261019.py`: label-only rendering from
  the saved extended statistics, without refitting or bootstrap recomputation.
  `SRC` is assigned and unused. The script does not reopen the example npz
  files and does not check the hash of `statistics.json`.

## Launch contract

These five files are frozen talk renderers for the home-directory package
`/global/homes/d/dkololgi/fmsu_tokyo_20261010`. They have no CLI. They do not
read `SLURM_JOB_ID`. Point them at the recorded inputs. Substituting a
corrected-label catalogue would rewrite the preserved slides.

| Script | Runs when imported | Output directory | Exclusive write |
| --- | --- | --- | --- |
| `desi_filament_visual_audit_20261017.py` | No. `main()` is guarded. | `…/filament_and_local_geometry_20261017` | No. `exist_ok=True`; JSON, npz, and figures are replaced. |
| `desi_slab_sensitivity_20261017.py` | Yes. No `__main__` guard. | same `data/` directory | No. Replaces `slab_sensitivity.json` and `slab_sensitivity_arrays.npz`. |
| `tokyo_wall_and_ssfr_overlay_20261017.py` | Yes | `…/talk_wall_overlay_20261017` | No |
| `tokyo_wall_and_ssfr_overlay_extended_20261017.py` | Yes | `…/talk_wall_overlay_extended_20261017` | No |
| `tokyo_galaxy_properties_labels_20261019.py` | Yes | `…/notation_corrections_20261019/F_galaxy_properties` | No |

Every script except the slab sensitivity pins `OPENBLAS_NUM_THREADS`,
`OMP_NUM_THREADS`, and `MKL_NUM_THREADS` to `1` at import. The slab script
pins only `OPENBLAS_NUM_THREADS`. The two wall overlays and the label
renderer also create their output directory at import, before the input
checks, so a later failure still leaves that directory behind. The slab
script does not create directories; it expects the audit `data/` tree.
`mkdir(exist_ok=True)` plus `savefig`, `write_text`, and `np.savez_compressed`
replace existing slides and receipts.

`desi_filament_visual_audit_20261017.py` does not write on import. Calling
`main()` re-reads `a06_desi/manifest.json` and `plotted_galaxies.npz`, checks
the npz hash, and re-reads those catalogue rows with `fitsio`. The shape statistic is
the 2D covariance contrast `(λ_major−λ_minor)/trace` of other displayed points
inside `R = 7/0.6766` Mpc, requiring at least 8 neighbours. The cut is the
80th percentile among interior valid rows, with the slab edge inset by `R`.
That cut is projected and descriptive. It is not the later 3D L/P/S statistic
and it is not a T-web class. PDF export requires the file to contain
`/FontFile2`. The audit JSON is written twice: once before the figures and
again after the crop counts are added.

`desi_slab_sensitivity_20261017.py` requires that audit npz and the frozen
`top20pct_cut`. It rebuilds neighbour trees inside `|height| < depth/2` for
depths 5, 10, and 20 Mpc. Fixed targets stay the original shape-valid rows
with `|height| < 2.5` Mpc. It does not redraw figures.

Both wall-overlay scripts hash-check each example npz against
`docs/evidence/filament_mock_examples_20261017/EXAMPLES.json` and hash-check
`density_checks_20261014/matched_galaxies.npz` against its `MANIFEST.json`.
The sphere test is `neighbours + centre`. `true_wall_sheet` is drawn in
`Δx`–`Δheight`; the string and the filament-error example are drawn in
`Δx`–`Δy`. The slide line `R = 7 h⁻¹ Mpc` is this same `7/0.6766` Mpc.
Diamond annotations are `Void`, `Wall`, `Filament`, `Cluster` at `P_WEB`
argmax. The colours are read from the theme keys `void`, `sheet`, `filament`,
`knot`, so the Wall diamond uses the sheet colour. Bar tick labels come from
`a07_combined/plotted_statistics.json` and are not hardcoded here.

The sSFR panel keeps a mass/redshift cell only when every trace bin and every
class has at least 20 rows in that cell. Both scripts draw 128 multinomial
resamples of sky regions with `default_rng(202610171)` and equal region
probability, then keep the 16–84% interval of draws that are jointly valid
for the binned curve and the class diamonds. The original script uses eight
quantile bins and requires at least 120 jointly valid draws. The extended
script appends a bin edge at trace `1.75`, sorts it into those quantile
edges, and requires at least 100 jointly valid draws. Its receipt describes
that extra edge as a split of the upper bin. `retained_n` and `common_cells`
are written at runtime; the recorded receipts are 63,347 galaxies in 74 cells
(original) and 59,774 galaxies in 69 cells (extended).

`tokyo_galaxy_properties_labels_20261019.py` only replots the extended
`statistics.json`, the original fraction panel, and `theme.json`. The trace
axis label is the expanded form
`Tr(T)=λ₁+λ₂+λ₃ ∝ δ`. Running it replaces
`A07_continuous_and_classes.{png,pdf,svg}` and copies the script and
`statistics.json` into the notation-corrections folder.

## Repository cleanup verification

All five newly tracked Python files passed syntax parsing; numerical pipelines
were not rerun. Three local `loa_page_*.log` run logs are ignored only under
`workflows/presentation/`; they remain on disk. The unrelated CONTEXT.md
label-provenance notice is left untouched and uncommitted.
