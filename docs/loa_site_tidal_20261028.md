# Loa Site tidal/clarity update

## Resumed release: historical displays retired at user request

The user approved a further one-node/32-CPU/desi-interactive/Scratch 20-minute
allocation, then explicitly instructed deletion of old displays. Allocation
60422798/nid200053 performs this bounded export/package/QA pass.

Ten historical displays were removed from the deployed asset tree and selector.
Their recoverable assets are at
`/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_site_retired_displays_20261028`;
scientific source products and prior Git/source versions are unchanged.
Current six-fit NGC/SGC density views and the local tidal page remain.

Lossless address/density coding was independently decoded in JavaScript against
every original integer, including the exported means. Exact comparison passed
for all twelve cap/fits and five displayed fields each. Mean residual coding
preserves the original quantized mean, not an approximate replacement.
Final browser QA passed all tidal controls, 390px layout, both density caps,
opacity and depth clipping, removal of historical options, and no browser errors
or external requests. QA directory suffix is `loa_site_tidal_20261028_qa_repack`.
Deployment assets are now 231,768,361 bytes before final archive compression.
Private deployment succeeded at https://loa-cosmic-web.astrodak99.chatgpt.site.
Tidal page: `/tidal.html`. Audience remains owner-only.
Final archive: 228,175,335 bytes (217.6 MiB), below the 256 MiB limit.
Pushed source: `d88b5b5375ac7201fc62b374eef2bee6f83beb0c`.
Deployment: `appgdep_6ae178870a148191b44d84dd1afd0bbc`.
Saved version:
`appgprj_6abdae1cd14c8191bf07f24c6b0a6a59~appgver_a102a70e89988191bb5fb6fee41014d0`.
Allocation was explicitly released after packaging and before deployment.

## Earlier blocked attempt (historical)

User requests: plot tidal fields, update the existing website, improve density
cloud clarity, and remove the survey-support display cutoff. Preserve actual
missing predictions rather than fabricate data. User approved one additional
CPU node, 32 CPUs, desi/interactive, Scratch, at most30minutes with release.

Site identity: appgprj_6abdae1cd14c8191bf07f24c6b0a6a59.
Existing audience verified owner-only custom access; no access change.
Checkout: /global/u2/d/dkololgi/loa_site_tidal_20261028_v8dHI9.
Final pushed source: 8ab2bc0b4da0d3157261c8e6fe8f788335e6e3fb.
Original single-commit content is preserved as
refs/backup/loa-tidal-unsplit-20261028; final trees are identical.

## Prepared and verified

- Current six fits in both caps re-exported from SHA-verified source HDF5,
  retaining every finite prediction-covered cell, no support cutoff. Existing
  2x arithmetic pooling and16-bit log-density precision retained (<0.0085%
  max relative error). No changed source fields or new inference.
- New tidal.html: all32768 cells in the216.512Mpc/h central NGC U-Net17 region,
  four draws; density/trace/ordered-eigenvalue slices and3D class boundaries;
  source widths649.5/974.3/1299.1Mpc/h and both exterior/footprint assumptions.
  Float32 export preserves threshold0.2 classes over all cells. These are
  provisional finite-context fields, not a validated full-survey tidal release.
- Stronger lighting, opacity slider, z-slab clipping, first draw rather than
  arithmetic mean as initial selection; galaxies available but initially hidden.
- Existing historical views retained with their original provenance/masks.
- Browser checks pass four draws, quantities, context/class controls, desktop
  and390px tidal layouts, both cap loads, opacity and clipping. No page errors
  or external requests. Screenshot/QA directory:
  /pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_site_tidal_20261028_qa_v1/.
  Select-option test ambiguity was resolved with explicit DOM value checks and
  change events. Mesh scene updates are serialized and reused.

Allocation60419859/nid200041 released after packaging, within30minutes.
Failed browser probes and packaging input attempts retained in task output.

## Publication blockers and safe continuation

First source push rejected at261.5MiB against64MiB per-push transport limit.
Recovered with eight <40MiB-content incremental commits, each fast-forward
pushed; final tree equals the tested original commit. No forced remote update.
Recovery helper: /global/u2/d/dkololgi/loa_split_site_push_20261028.mjs
(no credential stored; credentials entered hidden on stdin and used in memory).

Official source helper packaged:
/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_site_tidal_20261028_deploy.tar.gz.
Native private save/deploy rejected the398MiB archive against256MiB input limit.
No version/deployment ID was returned; live website remains unchanged.
Do not claim deployed. Do not retry this unchanged oversized archive.

Candidate lossless size reductions: derive arithmetic mean in browser from four
draws instead of shipping a fifth field; delta-encode monotonically increasing
cell indices before gzip; retain galaxy-only historical views while removing
their unused density payloads from the deployed assets, preserving original
products/history. Must revalidate masks, mean error, bounds and all views, and
measure archive size. Do not lower precision, hide cells or silently drop views.
Additional compute requires fresh resource-specific approval; the previous
allocation has been released. Refresh the source credential before more pushes.
