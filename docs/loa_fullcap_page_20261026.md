# Loa full cap Page update

Updated the private [Loa Page](https://chatgpt.com/space/page_ac5434c5e9e88191967a300a10ade530) on 26 October 2026. All six frozen fits are represented across 1,626 NGC and 846 SGC tiles. The Page is explicitly an independent-tile mosaic, not a joint reconstruction. Original pilot embeds and their saved interaction state are unchanged; historical headings now default to collapsed. The separate public Site and population Page were not modified.

## Coverage and display contract

- NGC: 3,293,489 retained galaxies, 3,010,154 property-eligible; 26 empty retained-galaxy tiles.
- SGC: 1,334,703 retained galaxies, 1,228,127 property-eligible; 23 empty retained-galaxy tiles.
- Total: 2,472 tiles, 4,628,192 retained galaxies; four draws per fit per tile.
- The density display uses native central z plane 8, four-draw arithmetic means and 4 by 4 in-plane pooling: 27.064 Mpc/h pixels, 6.766 Mpc/h thickness. It is not a full volume export. All contributing native support cells must be at least 0.5. Valid display pixels: 36,277 NGC and 15,748 SGC.
- Density colour quantization uses 16 log2 levels over [-4,4], with measured maximum relative error 20.303%. Origins are rounded to 0.01 Mpc/h; median redshifts to four decimals for display, while shell filtering retains original integer assignments.
- Per-tile class counts and six-model aggregate fractions are exact saved counts, not quantized density-derived labels. Counts describe retained galaxy-draws, not volumes or independent observations.
- No individual galaxy IDs or properties are exported. No tile, fit or outlier was removed.

## Sources and validation

The authoritative artifacts are under
/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1.

Output folder: page_fullcap_20261026_v2. PACKAGING.json records source hashes and each embed checksum; BROWSER_QA.json records passed local browser checks. Every archive checksum and consumed member checksum was verified. All tile identities and retained counts reconcile; all six per-cap class-count sums reconcile with SCIENTIFIC_DIAGNOSTICS_20261020.json.

Local browser checks exercised both caps, six fits, four redshift shells, Cartesian layer selection, tile detail, yaw and collaborator-state updates at widths 1200, 390 and 320. Screenshots were inspected at wide and narrow sizes. Optional preview-shell icon/tooltip CDN requests were blocked; no external network requests completed. Uploaded Page references have the expected byte sizes; in-Page rendering was not directly inspectable.

The native scientific review includes persistent transformer structure suppression, selection versus support, density outlier sensitivity, the existing four-versus-32 prefix comparison and the finite-domain mock limitations. Full-survey 32-draw expansion remains on hold. The CIGALE property explorer remains explicitly pilot-only.

## Execution and retained failures

User approved one CPU node, 64 requested logical CPUs, desi/interactive, Scratch, at most one hour. Allocation 60289725 on nid200043 completed and released after 7m32s (0.12556 CPU-node-hours); no GPU or inference was used.

The initial packaging attempt stopped at the 256 KiB embed limit (NGC 264,032 bytes); its empty output directory is retained. Reducing only displayed coordinate precision produced NGC 260,895 bytes and SGC 137,795 bytes. The first browser assertion flagged blocked optional preview-shell requests; all functional checks had passed. The corrected offline check passed. An attempted wrapper rerun refused to overwrite an existing generated preview; the final test reused the existing previews. The wrapper now explicitly allows replacement of these task-generated previews.

Required graphify refresh was attempted on the compute node and failed with filesystem flock OSError 524. No lock bypass or global-index refresh was performed. This source-index limitation does not change the passing scientific and browser checks.

## Publication

Page text update sequence 37; NGC embed 38; SGC embed 39; layout and historical disclosure 40. Original three embed references and their state were verified unchanged. No commit or push was requested.
