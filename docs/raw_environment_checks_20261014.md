# New raw environment measurements for the Loa sky wedge

2026-10-14. Requested scope: the same sky wedge as the previous plots, with all surrounding Loa neighbours retained. All requested statistics were calculated before applying CIGALE eligibility. No new ML inference was performed.

## Outcome

Independent density estimators corroborate the inferred property trends. Adaptive projected fifth-neighbour density has a stronger quenched-fraction contrast than the soft web classes, while 4- and 7-h^-1-Mpc counts have closely comparable amplitudes. Small fixed apertures are sparse: 76.1% of all supported query galaxies have no neighbour within 1 h^-1 Mpc.

Compactness retains a property association after adding large-scale density control. The quadrupole results depend on three-dimensional versus sky-plane geometry and have substantially smaller support. The local neighbour quadrupole and inferred gravitational shear proxy are distinct observables; their marginal property trends need not have the same sign.

## Measurements

Measured 103,788 supported galaxies at RA120-160 deg, Dec14.5-30.6 deg and 0.20<=z<0.30, using 352,198 Loa tracers in a padded Cartesian bounding box. Neighbours outside the display wedge are included. Query selection uses no CIGALE, density or inferred-class cut. All 5,436,413 Loa VAC rows were considered for the neighbour pool. Coordinates use Planck18 comoving distances; h=0.6766.

- **Projected Sigma5:** 5/(pi d5^2), where d5 is the projected comoving fifth-neighbour separation in h^-1 Mpc at the target distance. The window is |c dz/(1+z_target)|<=1000 km/s. Search radius is capped at10 h^-1Mpc;529 galaxies lack a fifth neighbour inside it and receive NaN plus an invalid flag. They are not assigned artificially low finite densities.
- **Spherical counts:** self-excluded redshift-space counts at1,2,4,7 h^-1Mpc. TARGETID exclusion is exact.
- **Compactness:** density inside2 h^-1Mpc divided by density in the2-7 h^-1Mpc shell: N2*(7^3-2^3)/[(N7-N2)*2^3]. No pseudocount. Undefined empty-shell values are flagged.
- **3D quadrupole:** M=mean(u u^T) for unit neighbour directions inside7 h^-1Mpc; A=(3 tr(M^2)-1)/2. The reported Q3=(N A-1)/(N-1) subtracts the independent-isotropic finite-N expectation. Negative estimates are retained.
- **Sky-plane quadrupole:** the same sphere neighbours are projected onto the target tangent plane; A2=|mean(exp(2i phi))|^2 and Qsky=(N A2-1)/(N-1). This removes the direct line-of-sight direction from the angular statistic but retains redshift-space neighbour selection.
- **Additional diagnostic:** LOS_P2=mean[P2(u dot line_of_sight)] is stored for future redshift-distortion checks.

Shape plots require at least20 usable directions for each quadrupole. The finite-N correction is exact under independent isotropic directions; it does not remove correlated clustering, mask incompleteness or redshift-space effects. A quadrupole detects one angular moment, not every possible form of anisotropy.

## Boundary handling and selection

The stored VAC boundary distance was calculated at a5-Mpc grid voxel centre. The new measurement instead looks up the frozen nside256 angular-support distance at each actual galaxy position and subtracts twice the maximum pixel radius. This is a conservative buffer for pixel-centre and boundary-pixel uncertainty. The sphere must fit within the angular support and the tracer radial limits; Sigma5 requires its measured d5 aperture and full velocity window to fit. This controls survey geometry at the mask resolution, not spatial variations in targeting completeness.

There are97,840 sphere-valid queries and99,502 Sigma5-valid queries. Requiring both plus current unique matched GALAXY CIGALE values gives96,258 rows with valid positive mass/SFR/AGE, nonnegative AV and FAGN in[0,1], with8<=logM<13. These are freshly selected current-CIGALE galaxies in the same geometric wedge; the historical property-cache ID restriction is not inherited. AGE denotes the supplied SED-fit parameter, and FAGN is a fitted contribution rather than AGN incidence.

## Controlled comparisons

Every curve within a suite/property uses identical galaxies, identical strata and identical population reference weights. Controls are0.1dex stellar mass and0.01redshift; the mass property controls redshift only. All four groups of every estimator require Kish effective N>=20 in every retained cell. Uncertainties are16-84% intervals from128 paired nside8 sky-block draws, retaining valid replicates; geometry has126 valid draws, density and density-controlled geometry128. Reference weights, group thresholds and cells are fixed during resampling.

Non-neural statistics and inferred trace/shear use four equal-total-weight quantile groups within each redshift bin. Tied values spanning a boundary receive identical fractional membership for every tied galaxy. This preserves equal bin weights without assigning arbitrary ranks to equal neighbour counts. Soft web groups retain P_WEB weights in void/wall/filament/cluster order. These different group definitions describe distinct environmental selections.

The density comparison retains91,066 galaxies in134 common mass/redshift cells. Stellar mass uses96,258 galaxies.

| Estimator | Quenched Q4-Q1 or cluster-void [pp] | 68% interval [pp] | Mean log sSFR contrast [dex] |
|---|---:|---:|---:|
| Projected Sigma5 | +9.518 | +9.153 to +9.783 | -0.2520 |
| 1 h^-1 Mpc counts | +3.296 | +2.835 to +3.642 | -0.0810 |
| 2 h^-1 Mpc counts | +6.463 | +5.963 to +6.863 | -0.1648 |
| 4 h^-1 Mpc counts | +7.819 | +7.387 to +8.150 | -0.2109 |
| 7 h^-1 Mpc counts | +7.679 | +7.370 to +7.885 | -0.2093 |
| Inferred trace | +6.346 | +6.016 to +6.590 | -0.1722 |
| Soft web classes | +7.420 | +6.976 to +7.677 | -0.1924 |

The1-h^-1-Mpc result should be read with its76.1% zero-count fraction: it is a very low-resolution environment indicator in this tracer sample. The median neighbour counts at1/2/4/7 h^-1Mpc are0/1/4/15. These sparse fixed apertures do not constitute a contradiction of stronger local-environment effects measured with adaptive projected neighbours.

## Geometry and the extra density control

Geometry requires >=20 neighbours and finite compactness;34,537 galaxies meet the property/geometry selection. Common mass/redshift support retains18,501 galaxies in47 cells. The additional large-scale-density control stratifies log10(N7) into0.2dex bins alongside mass/redshift, retaining4,539 galaxies in25 cells. At fixed narrow z, N7 is the directly measured large-scale tracer-density proxy. This remains a finite-width control rather than exact equality of density.

To distinguish the effect of adding this control from the change in sample, a separate mass/redshift-only baseline uses those exact4,539 galaxies and the original quantile memberships.

| Estimator | Mass/z:18,501 galaxies, quenched contrast [pp] | Mass/z: same4,539 galaxies [pp] | Mass/z/N7: same4,539 galaxies [pp] | 68% interval with N7 [pp] |
|---|---:|---:|---:|---:|
| Inner/outer density ratio | +7.271 | +5.997 | +6.233 | +3.851 to +8.140 |
| 3D direction quadrupole | -4.615 | -3.000 | -2.668 | -5.231 to -0.310 |
| Sky-plane quadrupole | -1.063 | +0.605 | +0.763 | -2.198 to +3.460 |
| Inferred shear proxy | +1.853 | +0.919 | +0.911 | -0.940 to +3.032 |
| Soft web classes | +5.953 | +4.559 | +2.168 | +0.439 to +4.642 |

Compactness is the clearest supplementary measurement: its roughly6-point quenched contrast is retained on the same galaxies after the density control. The3D quadrupole contrast remains negative but uncertain; the sky-plane quadrupole and inferred shear-proxy intervals include zero. Much of the difference between the full geometry sample and the density-controlled sample comes from the reduced support; do not attribute that entire change to density adjustment.

The same-row mass/redshift-only baseline has127 valid bootstrap draws; before/after interval differences are not a paired significance test for the change in contrast. No significance claim about that change is made.

The inferred shear proxy is q=sqrt(0.5*sum_{i<j}(lambda_i-lambda_j)^2) evaluated at posterior-mean eigenvalues. It is not the posterior mean of q. The 3D and sky-plane neighbour-quadrupole ranks have Spearman0.586; Q3 has Spearman0.380 with squared LOS_P2. The latter is partly a geometric relationship between components, not proof that redshift distortions cause the property trends.

## Files

- local_density_comparison: Sigma5 and1/2-Mpc/h counts versus soft web classes.
- multiscale_density_comparison:2/4/7-Mpc/h counts, inferred trace and soft web classes.
- geometry_comparison: compactness,3D/sky-plane quadrupoles, inferred shear and soft web classes.
- geometry_fixed_density: same geometry comparison with the additional N7 control.
- Eight *_properties figures: each measured statistic against all seven properties on its physical x-axis.
- raw_environment_catalogue.npz: per-galaxy quantities, TARGETIDs, coordinates, CIGALE fields and validity masks.
- statistics.json / contrasts.csv: all seven properties, means, intervals and correlation diagnostics.
- analysis_membership.npz: exact retained rows per suite/property; indices refer to raw_environment_catalogue.npz.
- measurement_provenance.json, analysis_provenance.json, additional_QA.json, tie_membership_audit.json, MANIFEST.json and scripts/: reproducibility and validation.

## Verification

The bounded neighbour measurement took13.9seconds on one thread after the initial FITS-byte-order fix. User requested login-node execution; no Slurm job was submitted, cancelled or modified.

Passed synthetic exact-self-exclusion, known fifth-neighbour distance, rotation invariance, and isotropic finite-N tests. Passed tied-value equivalence, equal group-weight, preserved-signal and constant-estimator-null tests. Twelve deterministic real galaxies were checked by brute force against all padded tracers, including projected candidate completeness and fifth-neighbour distance. Counts are monotonic with radius; identities are unique; quadrupole bounds and valid Sigma5 finiteness pass.

The recomputed R7 counts exactly match all94,415 overlapping historical cached galaxies. Support maps were freshly SHA256-verified against the frozen input receipt. Full FITS checksums are inherited from frozen receipts, not recomputed; source row count and exact current-CIGALE TARGETID/redshift parity were checked. Derived arrays and deliverables have fresh hashes.

No neural inference, training, posterior resampling, release overwrite, CFS mutation, Git commit or push. graphify is unavailable; targeted source inspection was used. Scientific figures retain the established dark presentation style and are available as PNG/PDF/SVG.

Intervals describe spatial sampling only. Common CIGALE parameter errors, variable angular completeness and redshift distortions are not removed by using non-neural environment statistics. The new plots supply an independent-estimator consistency check within the selected wedge.
