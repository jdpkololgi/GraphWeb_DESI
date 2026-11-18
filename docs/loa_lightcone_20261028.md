# Full-redshift Loa density wedges

Subsequent extension: [full-cap tidal wedge display](loa_lightcone_tidal_20261028.md)
adds threshold-controlled T-web layers to the companion site. The density-only
contract below records the original version; density and galaxy data are unchanged.

User requested continuous structures across the available redshift range, rather
than isolated cubes, with an observed BGS_BRIGHT overlay. Website page:
`lightcone.html`, in existing owner-private Site
`appgprj_6abdae1cd14c8191bf07f24c6b0a6a59`.

## Contract

- Existing shared-state expanded U-Net seed17 fields, NGC and SGC, four saved
  draws; no training, new sampling, truth reads or modifications to inference.
- Full available field range 0.15 <= z < 0.55. Radius is comoving Mpc/h, with
  redshift ticks. Distance mapping is recovered from the actual validated parent
  coordinates; maximum interpolation error 1.36e-5 Mpc/h. Native cells 6.766
  Mpc/h; learned density target Gaussian R7. No extra smoothing or fine-resolution
  scientific claim.
- Native-cell lookup on half-degree RA and <=6.766 Mpc/h radial intervals, at
  one-degree declination-bin midpoints. Selected declination bins are averaged
  with solid-angle weights. This is angular quadrature, not a Cartesian slab.
  Finite angular sampling does not establish exact coverage between samples.
- A field sample requires prediction coverage and positive pooled random
  support. Rays missing any selected declination-bin sample stay blank. No
  interpolation across masked gaps; an alternate layer shows the fraction of
  sampled angular weight with support. This finite-resolution mask is not an
  exact survey boundary or a calibrated reliability probability.
- Full observed Loa v2.1 `BGS_BRIGHT_full_HPmapcut.dat.fits` selection:
  `(BGS_TARGET & 2) != 0`, ZWARN=0, DELTACHI2>=25, GALAXY, and stated redshift
  interval. Unique TARGETIDs; no thinning, magnitude recut or field-support cut.
  NGC 3,868,333 and SGC 1,568,080; total 5,436,413. Galaxies are selected by the
  same declination edges and RA span. Points have no inferred T-web class.
- Density only: local tidal boxes are not stitched into full-cone T-web labels.
  Apparent connectivity is exploratory, can result from projection overlap, and
  is not calibrated filament recovery or a connectivity statistic.

## Implementation and checks

Exporter: `workflows/presentation/loa_lightcone_20261028.py`.
Output: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_lightcone_20261028_v1`.
`manifest.json` binds field/parent/exporter hashes, catalogue size/mtime,
per-slice TARGETID hashes and compressed asset hashes. Catalogue provenance is
size/mtime, not a full-file hash; before/after source identity is checked.

`QA.json` passes source hashes, positive finite densities, exact byte roundtrips,
coordinate-direction/cap agreement, monotone distance mapping, unique galaxy
IDs and exact partition totals. 91 NGC and 55 SGC one-degree bins have field or
galaxy data; unpopulated declination bins remain empty.

Site `test-lightcone.cjs` passes synthetic wraparound, half-open RA selection,
solid-angle weighting, draw averaging and missing-support tests.
`browser/QA.json` passes real-data native reference pixels, four-draw means,
exact full-RA galaxy counts, both caps, slice/draw/layer/overlay controls,
390px layout and PNG export, with no page errors or external requests.
Desktop image inspected: full range and density/galaxy overlay render visibly.

Approved allocation 60425806/nid200007: one CPU node, 32 CPUs, desi/interactive,
Scratch, up to 30 minutes, release on completion. No GPU or new inference.
Publication and allocation cleanup are recorded after completion, not inferred
from local test success.

Required AST graph update attempted on allocation step 3; it failed before
extraction with the existing `fcntl.flock` OSError 524. No lock bypass or global
graph refresh from stale output was performed.

Allocation released after 9m13s: allocation COMPLETED 0:0; schema, export and
browser steps COMPLETED 0:0; graph refresh step FAILED 1:0 as above. Native
host launch was explicitly reviewed because this saved chat's old sandbox
failed before Slurm invocation. No chat/security configuration was changed.

Site source upload exceeded the 64 MiB per-push transport limit on first try.
Recovery splits identical final content into <40 MiB commits and preserves the
original local commit under a backup ref. No galaxy thinning or field reduction.

The combined archive then exceeded the separate 256 MiB deployment limit.
The existing site's full-cap/cube/tidal assets are preserved. Light-cone assets
were moved intact to companion owner-private Site
`appgprj_6ae18b46013c819183a767bbd43887f2`, with its own root page and navigation
back to the existing viewer. The original `lightcone.html` route redirects to
the companion. This changes packaging only, not the data or scientific scope.

## Published result

Both native deployment calls succeeded on 2026-10-28:

- Companion: https://loa-lightcone-slices.astrodak99.chatgpt.site
  source `340e4bd8f61f0a6a359bc8dbd3b2a844c7b56df5`, deployment
  `appgdep_6ae18c5243f08191845dae2a5696ba2e`.
- Existing site's navigation/redirect: https://loa-cosmic-web.astrodak99.chatgpt.site
  source `4ca330fc5eb20d4be601bd3b492d472e4da98ed3`, deployment
  `appgdep_6ae18c8d06e48191a460679acc5b91b2`.

Owner-only access preserved. Companion source is
`/global/u2/d/dkololgi/loa_lightcone_site_20261028`. Renderer and math scripts
were compared byte-for-byte after relocation; generated data moved intact.
Browser evidence precedes relocation and a small legend-label clarification;
hosted deployment success is verified, but a separate authenticated production
browser test was not run. Original site's historical test script stages into its
old `dist/lightcone-data`; use the companion root for future light-cone tests.
