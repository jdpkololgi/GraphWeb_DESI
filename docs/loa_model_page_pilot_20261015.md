# Loa field model and galaxy property Page pilot

## Scope and status

Observed-Loa inference is complete for two finite domains, three field-model
families and two seeds, four draws per fit/domain (48 draws). All twelve result
receipts pass numerical checks. Browser verification passed and all three
interactive embeds are saved in the private Page. This is neither a full-survey inference campaign nor calibration
qualification. No training, mock truth, oracle coarse field or automatic
inference retry was used.

Private Page: <https://chatgpt.com/space/page_ac5434c5e9e88191967a300a10ade530>.
The original Loa Cosmic Web Site is unchanged.

## Selection and comparability

The parent atlas plan contains 1,626 north and 846 south tiles, with
5,436,413 VAC rows across supported redshift 0.15–0.55. The pilot selects the
largest observed-row-count tile in each cap, before examining galaxy properties
or inferred outcomes. It is a throughput/property-rich pilot, not a
representative survey sample.

| Cap | Tile | Retained galaxies | Actual galaxy redshift |
| --- | --- | ---: | --- |
| North | NGC_006_010_005 | 26,982 | 0.150030–0.215687 |
| South | SGC_001_014_004 | 23,808 | 0.150001–0.221140 |

Half-cell interpolation edges are excluded. Expanded observed channels reproduce
the legacy central observations exactly in both cases. Catalogue joins require
unique TARGETID and exact RA/Dec/redshift agreement with the brightness source.
The parent redshift range must not be described as the coverage of these cores.

The earlier CFM, expanded U-Net clean endpoint and transformer velocity each use
seeds 17 and 29. The four expanded fits use frozen step-9984 EMA checkpoints;
legacy coarse/fine fits use their existing frozen checkpoint pair. All models
share the same retained physical 32×16×16 cores. Coarse draws are shared across
families within each seed. The expanded models share fine noise; the legacy
fine domain has a different shape and its own addressed noise.

The coarse and legacy samplers use 64 Heun steps (128 network evaluations);
the expanded sampler uses NFE128. Tides use a common finite free-space
local-plus-wide operator, with mean density beyond the wide domain. This is
not the learned full-cap exterior used by a different existing visualization.
Saved eigenvalues are derived from each model's own density/tidal draw.

## Numerical and property checks

Positive density, finite outputs, coarse-mass conservation, free-space core
trace, interpolated galaxy trace, source bindings and output checksums pass.
Maximum relative coarse-mass error is below 6e-16 and maximum core trace error
below 9e-15. These are engineering checks, not evidence of calibration.

For the property explorer, require unique unambiguous CIGALE GALAXY matches,
saved catalogue support plus interpolated core support ≥0.5,
9≤log10(Mstar/Msun)<12 and corrected r<19.5. North photometry receives the
registered −0.04 mag correction. This leaves 25,425 north and 22,943 south
galaxies (48,368 total). All seven displayed properties are valid for that
selected population in this pilot; the implementation retains property-specific
validity rather than assuming this for future data.

The explorer shows log stellar mass, log SFR, log sSFR, log stellar age in Gyr,
AV, AGN fraction and quenched fraction (log sSFR<−11). Invalid/nonpositive SFR
does not silently count as quenched. For λ=0.2, use the fraction of field draws
in each environment as weights. P12-A is a separate galaxy-level comparator,
not a field or a matched-exterior posterior.

Stored aggregates are Σw, Σw² and Σwy by cap/model/mass/redshift/brightness/
environment/property. Means and effective counts are recomputed after filters.
All 98 model/cap/property weight sums reconcile to the eligible property counts.
Effective N is a weight-concentration diagnostic, not independent information or
an uncertainty interval. No causal environmental interpretation is claimed.
Only aggregates and a deterministic unclassified coordinate sample are embedded;
individual TARGETIDs and CIGALE measurements are not published.

## Presentation and validation

Source templates are under `workflows/presentation/loa_model_page_20261015/`.
Each cap gets a separate self-contained HTML embed to respect the 256 KiB Page
limit. Within each embed, model views share rotation, physical bounds and
isosurface thresholds. Density draws are pooled by two to 16×8×8 display cells
(13.532 Mpc/h), with <0.009% log-quantization error. Environment fractions
retain exact 1/32 increments after four-draw/native-cell pooling. Entropy is
pooled native-cell posterior entropy, not entropy of the pooled fractions.

The property embed provides all seven properties, two simultaneous model
comparators, selection filters, a mass-stratified view, exact values and CSV.
Empty selections are not filled or extrapolated. External requests are forbidden
in the browser test. Source tests also compare Python-decoded aggregate means,
weights and effective counts against JavaScript output.

## Execution and retained attempts

- CPU preparation: job 59820724, completed, preparation 269.68 s; released.
- GPU inference: job 59827704, inference receipt passed in 295.60 s. Subsequent
  packaging exceeded the combined-embed size limit and the step exited 1;
  the completed inference is retained and was not repeated. Allocation released.
- CPU packaging: job 59827975, split embeds generated with 98 reconciliations;
  browser check caught a missing JavaScript arrow. Allocation released.
- Final approved CPU-only browser verification: job 59828019, ten-minute limit;
  all browser checks passed, wide/narrow screenshots inspected, then released.

Page readback verifies the updated title and three full-width visualization
blocks. The exact uploaded HTML was tested locally; the live Page sandbox was
not independently inspected through a browser tool. The final queue is empty
and no tmux sessions remain. `sacct` and `graphify update` remain blocked by
the session's bwrap temporary-directory error; no permission settings or
security hooks were changed. Do not confuse the failed packaging steps with
the separately passed inference receipt.

Scratch root:
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_model_pilot_20261015_v1`.
`PREPARED.json`, `INFERENCE_COMPLETE.json`, per-cap/model `COMPLETE.json` and
`draws.npz` are authoritative. Earlier packaging directories are retained,
not overwritten; the current corrected package is `page_v3/`.

## Scale-up estimate and decision boundary

Measured mean fine sampling plus tide time is about 7.26 s/draw for the older
CFM, 3.65 for U-Net and 1.18 for transformer; coarse sampling is about
5.77 s/draw/seed. The largest measured allocated tensor memory is 1.87 GiB.
These timings come from the two dense pilot domains on one A100, not a full
footprint benchmark.

Reusing coarse draws across all three families, a simple independent-tile
projection over 2,472 tiles gives 98.2 GPU-hours for all six fits at four draws;
projecting total pilot elapsed time gives 101.5 hours. Linear 32-draw projection
is approximately 785–812 GPU-hours. This excludes survey-wide preparation,
queue time, storage/transfer variability and any new global/exterior method.
It is an engineering scenario, not a production budget or global-coherence
claim. No full-survey allocation is authorized or launched.
