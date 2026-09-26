# Independent alignment controls while awaiting Jade

Owner: GraphWeb_DESI, W1/W2 of the Mock-to-LOA Alignment Investigation.
No mock regeneration, calibration changes, new phase exposure or VAC replacement.
All catalogue reads are read-only. Existing exposed ph006 supplies the original
observed-stage mock; producer velocity checks use already-exposed ph000 only.

## Boundary and tile-count controls: completed

Reproduced the earlier Loa and original observed-mock histograms exactly before
changing the footprint. Removed one and two neighbour rings from the original
common occupied NSIDE256 sky, using NEST ordering and consistent child pixels.
Missing HEALPix neighbours at exceptional vertices are ignored, not treated as
survey holes. The five focused tests include this case and redshift-bin migration.

| Common sky | Loa/mock, z=.45–.55 SGC | Loa/mock, z=.45–.55 NGC |
|---|---:|---:|
| Original | 1.8565 | 1.8395 |
| Remove one pixel ring | 1.8375 | 1.8501 |
| Remove two pixel rings | 1.8231 | 1.8555 |

The high-z deficit survives removal of boundary regions. Low-shell two-ring
ratios are0.9841/0.9865. The two-ring region retains34,909/82,291 high-shell Loa
and19,148/44,350 mock galaxies, so this is not a tiny tail sample. This argues
against the mismatch arising predominantly from these footprint-edge pixels.
It does not certify exact angular-mask parity or exclude interior selection
variations. NTILE strata are stored separately; equal NTILE alone does not imply
equal targeting density, fibre competition or completeness.

## Paired velocity-convention diagnostic: completed

First10,000 raw rows of each producer N/S and canonicalv1 file, CEN=1, selected
on .15<=Z_COSMO<.55. Matches use exact HALO_ID/RA/DEC and verify identical
Z_COSMO. This prefix sample has matched pairs only in .35–.55; it is not a
representative full-redshift census and does not establish galaxy-ID semantics.

A simple v_new=s(z)*v_old cannot account for all differences. Even allowing
an independent best scalar for every matched central leaves component RMS
residuals142–170 in the stored velocity units. Median vector cosines are
0.873–0.916. This rejects a purely multiplicative explanation for these pairs;
it does not identify the physical central-velocity or random-assignment recipe.

All three files numerically obey the usual radial relation with c approximately
300,000: (Z-Z_COSMO)/(1+Z_COSMO) ~= v_los/300000. Effective fitted c is within
0.035 of300000; residual RMS is0.0087–0.0089 in velocity units. Using299792.458
instead leaves0.276–0.316 RMS. This small shared convention does not explain the
large inter-version velocity differences or the count deficit. No corrections
were applied. These are numerical fingerprints, not executed-source provenance.

Evidence: `evidence/p12a_independent_controls_20260926/VELOCITY.json`,
`RSD_SPEED.json`, and `boundary/COMPLETE.json`.

## DR1 / Loa release and random-footprint control

Completed on CPU58921665/nid004148. Uses ironv1.5 versus Loav2.1 full BRIGHT,
ZWARN0, DCHI>=25, GALAXY, .15<=Z_not4clus<.55. No weights. Index0 randoms define
NSIDE256/512 shared support and a >=16-randoms-per-pixel interior sensitivity.
No assertion of exact masks or unbiased area from finite random occupancy.
Matched-ID gains/losses, redshift migration and catalogue-only IDs are separate.
Neither membership differences nor quality changes uniquely identify additional
observations versus changed reductions/masks; exposure histories are not joined.


Tile-stratified two-ring high-shell ratios (SGC/NGC): NTILE<=1=2.139/2.099,
NTILE2=1.912/1.944, NTILE>=3=1.664/1.740. Thus the discrepancy is not confined
to singly covered regions; different observation completeness within strata
remains possible. No covariance or independent-region significance was estimated.


## DR1 results and interpretation

6,405,170 full-catalogue TARGETIDs are shared. On random512 shared support,
2,685,453 galaxies pass science cuts in both releases; none has |delta r|>0.01mag,
and801 have |delta z|/(1+z_DR1)>.005. There are689,136 matched-ID gains, of which
686,381 (99.6002%) had LOCATION_ASSIGNED=false in DR1. Remaining exclusive
first-failure reasons are355 ZWARN,849 SPECTYPE,1,067 DELTACHI2 and484 redshift
range. This supports added assignment/observations as the main origin of the
matched-ID gains, rather than wholesale photometric changes or merely lowering
quality thresholds. It is conditional on these catalogues and shared sky, not a
redshift-success model for every unobserved target.

| Loa/DR1 counts | .15–.25 | .25–.35 | .35–.45 | .45–.55 |
|---|---:|---:|---:|---:|
| Random512 shared sky, SGC | 1.543 | 1.547 | 1.518 | 1.492 |
| Random512 shared sky, NGC | 1.357 | 1.350 | 1.331 | 1.314 |
| Random512 interior, SGC | 1.455 | 1.457 | 1.430 | 1.410 |
| Random512 interior, NGC | 1.299 | 1.291 | 1.274 | 1.258 |

The increase is slightly smaller at high redshift, unlike the sharply growing
Loa/original-mock ratio. A simple DR1-versus-Loa observed-density difference does
not reproduce that redshift dependence. This is not a matched Y1-mock validation.
A switch of primary data release is not indicated by these controls alone.

The original Loa/mock common mask covers sky outside DR1. Its Loa/DR1 ratios
~2 include added area, and must not be interpreted as density ratios on matched
DR1 sky. Random-supported/interior controls reduce this geometric contribution.
One finite random catalogue per release and pixels still leave subpixel mask
uncertainty; exact tile/veto-map replay is open. Increasing random resolution
and removing boundary pixels change normalizations, not the qualitative slope.

## Validation and reproducibility

Five focused unit tests pass. Every saved broad-shell count reproduces compact
fine-bin histograms; source, helper, mask and output hashes verified. Exclusive
ID partitions and gain categories conserve each cap/fine-z count; accepted
redshift migration is explicitly included in the gain/loss balance. The boundary
baseline exactly matches the prior observed-mock/Loa census. Both figures were
visually inspected, including common y-scales within comparison rows.

Figures: [release/boundary](figures/p12a_independent_controls_20260926/release_boundary_controls.png),
[tile strata](figures/p12a_independent_controls_20260926/tile_strata_controls.png).
Compact evidence and validation: `evidence/p12a_independent_controls_20260926/`.
Allocation released after completion. No model refitting or correction performed.

Next: obtain Jade's intended targeting/population prescription; complete exact
angular/observation recipe parity before running her parent through the Loa chain.
These controls deprioritize broad release mismatch and footprint-edge artefacts;
they do not establish a unique LF/K/evolution cause or qualify a replacement VAC.
