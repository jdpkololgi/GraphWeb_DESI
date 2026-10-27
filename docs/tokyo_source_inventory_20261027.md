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

## Repository cleanup verification

All five newly tracked Python files passed syntax parsing; numerical pipelines
were not rerun. Three local `loa_page_*.log` run logs are ignored only under
`workflows/presentation/`; they remain on disk. The unrelated CONTEXT.md
label-provenance notice is left untouched and uncommitted.
