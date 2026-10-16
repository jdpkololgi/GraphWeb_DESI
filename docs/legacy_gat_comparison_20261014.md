# Old GAT catalogue and weak property trends

The same-source CIGALE comparison is now complete. The old sample has stronger
property contrasts, but matching stellar-mass distributions and label treatment
substantially narrows the gap. This does not establish that a historical
mass-completeness correction caused the difference, or that GAT is more accurate.

## Established so far

Found `/global/homes/d/dkololgi/GraphWeb_DESI/DESI_BGS_PRERELEASE_VAC.pkl`:
199,314 rows, 198,827 distinct integer TARGETIDs, saved GAT labels and four
probabilities. Probabilities sum to one and their largest component reproduces
every saved label. Median largest probability is 0.9281. File modification date
is January 26, 2026; this does **not** establish that it is the exact catalogue
behind the September 2025 slide. The archived inference notebook identifies
TNG300 snapshot 99 and a 1e9 solar-mass simulation cut, but contains several
checkpoint variants. Do not assign a unique checkpoint without further evidence.

The old redshift range is 0.0100293–0.05999996. The frozen P12-A CIGALE full VAC
has 5,436,413 rows at 0.15000002–0.54999828. Seven integer TARGETIDs occur in
both, but **none** passes |delta z|/(1+z)<0.001: discrepancies are 0.186–0.331.
There is no usable same-galaxy comparison in these saved products.

The z=0.15 lower bound traces to a **training-label failure found July 13–14,
2026**, recorded in Illustris `SCIENCE_LOG.md` under “SHELL-0 FAILURE SOLVED”,
the subsequent snapshot-shell correction, and “Remediation P1+P2+P3 done”.
In the audited low-z mock wedge, 100% of z<0.15 rows had BOX_INDEX=-1, versus
0% above the boundary. Their labels were as unpredictable from position as
shuffled labels (reported R² 0.004 versus -0.002). The record attributes this
to invalid mapping of the low-z shell onto the z=0.2 T-Web box; the exact mock
stitching boundary was empirically diagnosed, not independently established
from a recovered upstream generation configuration. The July 14 valid-cache
rebuild excluded that shell and required BOX_INDEX>=0 for supervised rows.
Later full-footprint/P12-A work inherited science shells
[0.15,0.25,0.35,0.45,0.55]. Current preparation includes z=.1–.6 context but
only emits predictions for the four science shells. This is a historical
label-validity restriction, not a stellar-mass-completeness or physical limit.
Restoring overlap with the old VAC requires valid low-z mock labels and
training/validation there; lowering an inference cut is insufficient.

487 old IDs each have two rows; all have differing mass/redshift/probabilities,
and 53 have different hard labels. Reproduction plots retain the saved rows;
the controlled colour check and planned CIGALE join exclude all repeated IDs.

The saved pickle contains FastSpecFit LOGMSTAR, fluxes and absolute magnitudes;
it contains no CIGALE properties. The current legacy loader applies only
ZWARN==0, DELTACHI2>=25, LOGMSTAR>=9 and BGS_TARGET!=0. Applying those four cuts
to the surviving 768,700-row `loa-combined-lowz.fits` gives **exactly 199,314
rows and the same TARGETID multiplicities as the pickle**. The separate
mass-completeness notebooks exist, but an extra per-galaxy cut from them is not
evident in this saved VAC. **This does not disprove that the completeness
calculation was used to choose the mass threshold and z<=0.06 ceiling.** The
earlier version of this report incorrectly treated those possibilities as
exclusive and overstated what the exact sample match established.

The notebooks use a redshift-dependent mass limit: the 90th percentile of
limiting masses among the faintest 20%, separately for bright/faint targets.
The saved plot in cell 35 of `2_loaSelectedSample_massCompletenessCut_forDK.ipynb`
shows the bright-sample limit crossing log mass~9 around z~.05–.06 (visual
estimate, not a recovered numerical fit). Thus the user's recollection is
consistent with the surviving calculation: estimate the completeness curve,
choose a redshift ceiling near where it meets the target mass, then use a
fixed mass cut. Commit `8848f4d` (June 17, 2025) changes the loader's upper
redshift from .1 to .06 alongside adoption of Manasvee's column definitions;
it does not explicitly state the rationale. The notebook source and saved
outputs also reflect different selections, so they are not a clean run receipt.

There is a mass-unit qualification: the notebook converts FastSpecFit h=1
masses using h=.704; the loader cuts raw LOGMSTAR>=9. Under that conversion,
the loader threshold is log physical mass~9.305, about 2.02e9 solar masses.
This is not exactly the same physical 1e9 threshold as intended for TNG.
The saved bright completeness curve at z=.06 is approximately 9.1 in the
converted units, below that converted loader threshold. These approximate
plot readings support the selection rationale but do not certify completeness
of every saved row or establish the exact original decision. The original
numerical completeness output under Manasvee's CFS directory is permission-
denied; no access workaround was attempted. `completeness_notebook_saved.png`
is the unchanged embedded notebook image, extracted for inspection.

## North/south and prediction-row ordering audit

The exact historical bug fix is commit `5296648` (August 7, 2025): reorder
`zcat` to `[north, south]` after generating hemisphere coordinates. Fresh
reconstruction of the saved old selection reproduces all 199,314 pickle rows
in that corrected order, including TARGETID, RA, DEC, Z and LOGMSTAR exactly
(139,545 north, 59,769 south). None of its TARGETIDs matches the unreordered
input at the same row. This confirms the saved catalogue uses the corrected
order; it does not independently replay the old model.

For current P12-A, the frozen production sources match the inspected code.
Within each core the same input row indices select IDs, coordinates, positions
sampled from the U-Net, posterior contexts and output rows. Tensor coordinate
conversion changes spatial-axis order, not galaxy order. Whole structured
records are sorted by TARGETID during merging. CIGALE joins use IDs plus
position/redshift checks, and plotting caches preserve original VAC row indices.

Fresh bounded audit (`ordering_audit.json`):

- All 5,436,413 output IDs exactly match the selected input IDs, uniquely.
- Eight small cores span both caps and all four redshift shells: 1,009 saved
  rows have exact input IDs/sky coordinates and correct spatial ownership;
  sky-to-Cartesian reconstruction is exact. Saved posterior draws reproduce
  their probabilities and mean eigenvalues exactly; all original output
  fields agree with the merged VAC, accounting for its documented quality bit.
- 256 rows sampled across the current plotting cache agree exactly with the
  original and CIGALE-augmented VAC. All 255 matched rows reproduce five CG15
  properties from fresh source-row reads; sky separations are zero and maximum
  normalized redshift difference is .000253, below the .001 join threshold.

No evidence of the old ordering failure was found. The ID census is complete;
the fresh prediction/property checks are sampled, not a full model replay or
proof that no upstream scientific error exists. Existing full-run receipts
also record complete saved-draw/summary validation at the original merge.
Reproducer: `workflows/presentation/ordering_audit_20261014.py` (one CPU thread,
no inference, no full CIGALE scan).

## Reproduced plots and first numerical check

Outputs: `/global/homes/d/dkololgi/fmsu_tokyo_20261010/legacy_gat_comparison_20261014/`.

- `legacy_colour_mass.png/pdf`: observed colour, -2.5 log10(FLUX_G/FLUX_R).
- `legacy_rest_colour_mass.png/pdf`: ABSMAG01_SDSS_G minus ABSMAG01_SDSS_R,
  SDSS bands shifted to reference z=0.1. Its red peak near 0.94 resembles the
  supplied slide more closely than observed colour, whose red peak is ~0.83.
  This resemblance is suggestive, not proof of the slide's exact origin.
- `audit.json`: source hash, schema checks, IDs/redshift discrepancies, counts.
- `selection_check.json`: exact saved-sample check.
- `legacy_colour_control.json`: compact numerical check, without error bars.

The figures separate counts from separately normalized class distributions
and use one intensity scale across the four colour–mass panels. A strong
blue/red bimodality within every panel is not itself an environment effect;
the changing relative numbers in the two populations are the relevant signal.

For the 198,340 nonrepeated-ID rows, mean rest-frame cluster-minus-void colour
is +0.12542 mag. Balancing 0.1-dex mass / 0.01-redshift cells leaves +0.10392
mag, retaining 193,679 rows in 98 cells with at least 20 objects in every class.
The raw and balanced summaries have slightly different support. No sampling
uncertainty was estimated for this quick diagnostic. This is an old-catalogue
association, not a measured old-minus-P12-A difference.

## Completed same-source CIGALE comparison

The user explicitly approved one bounded login-node run after Slurm inspection
failed. Scanned the 33,346,076-row, 20 GB Zou source once, one CPU thread.
Of 198,340 old rows with nonrepeated IDs, 193,706 have a unique position/redshift-
consistent CIGALE observation. Another 4,634 IDs have multiple consistent
observations and are excluded. Fourteen source rows fail position/redshift checks.
The seven-property selection retains 191,951 old and 102,623 current P12-A
galaxies: valid CG15 mass, SFR, fitted age, attenuation and AGN fraction, galaxy
spectral type, and 9<=log10 mass<12. Current P12-A uses the supported
RA120–160/Dec14.5–30.6/z.20–.30 wedge, without the neighbour-edge cuts from
the preceding density investigation. Its probabilities and current CIGALE
measurements come from the SHA-verified `raw_environment_catalogue.npz` cache.

All seven property panels compare labels and probability weights on identical
rows and common mass/redshift cells **within each sample**. This retains
187,358 old rows and 75,321 P12-A rows for non-mass properties. Mass controls
redshift only and retains all property-valid rows. All other panels use
0.1-dex mass / 0.01-redshift cells with >=20 effective rows in each class for
both label and probability methods. The mass distribution still differs between
the two samples; balancing classes within a sample does not match samples.

A second check therefore matches the mass distribution **between samples**.
The shared 0.1-dex bins span 10.2<=log10(M/Msun)<11.3. Each bin's reference
weight is proportional to the smaller eligible count in the two samples.
Within each mass bin, each sample retains its own pooled eligible redshift
distribution, applied identically across all four environments and both
membership methods. This retains 63,373 old and 74,816 P12-A rows.

| Comparison | GAT quenched contrast | P12-A quenched contrast | GAT log-sSFR contrast | P12-A log-sSFR contrast |
|---|---:|---:|---:|---:|
| Each sample's mass mixture, probability weights | +21.87 pp | +8.97 pp | -0.518 dex | -0.227 dex |
| Shared mass mixture, probability weights | +14.37 pp | +8.81 pp | -0.326 dex | -0.222 dex |
| Shared mass mixture, hard labels | +14.89 pp | +11.85 pp | -0.337 dex | -0.294 dex |

All contrasts are cluster minus void. Quenched means log10(sSFR/yr^-1)<-11.
All stages balance redshift across environments separately within each sample.
The numerical gap shrinks from 12.90 to 3.04 percentage points across these
choices, but stages change support and weights: this is not a causal fraction
explained or a same-galaxy model test. The final hard-label intervals are
14.16–15.69 pp for GAT and 11.30–12.78 pp for P12-A (16–84% spatial resampling).
Intervals include sky-sampling variation only, not SED or transfer systematics.
The models have different confidence distributions: median maxP=.9281 old,
.7852 current on their property-valid samples. More confident does not mean
better calibrated; hard labels are a diagnostic, not a preferred truth estimate.

The user's suggestion about the mass distribution is therefore supported:
which stellar masses we average over matters materially. The saved completeness
curve supports the intended low-z/fixed-mass selection, although exact mass
conventions and the original numerical cut need reconciliation. The surviving
loader and saved IDs establish the implemented cuts, not why they were chosen.
The much lower redshift allows fainter, lower-mass galaxies into a brightness-
limited survey. A completeness test using consistent mass units remains to be
done. Matching mass also leaves a substantial
offset in the overall sSFR/quenched levels between samples; redshift, colour-
dependent selection at fixed mass, and SED behaviour need separating before
interpreting that offset as evolution or model performance.

New plots in the output folder:

- `cigale_properties_raw_same_rows.png/pdf`: seven properties, raw means on
  the rows retained for each sample's corresponding balanced comparison.
- `cigale_properties_controlled.png/pdf`: seven properties with mass/redshift
  balanced within each sample; mass panel balances redshift only.
- `cigale_mass_ssfr_distributions.png/pdf`: four classes for each catalogue,
  probability-weighted class fractions and one shared intensity scale.
- `selection_and_weighting_check.png/pdf`: the three comparisons above.
- `legacy_cigale_colour_mass.png/pdf`: current CIGALE mass and observed flux
  colour on valid old-sample rows; not an exact remake of the rest-frame slide.
- `legacy_cigale.npz`, `join.json`, `comparison.json`, `common_mass_check.json`,
  `validation.json`: extracted data, provenance and numerical results.

## Follow-up status

The subsequent [consolidated investigation](legacy_gat_sequential_tests_20261014.md)
contains the executed selection, confidence, CG5, label-scale, checkpoint,
density-thinning and low-z provenance tests, plus the simulation alternatives
and pretraining design. It supersedes the pending status in the list below.
In particular, eight-query replay now identifies the January 15, 2026 GAT
checkpoint. Public source shows separate faint low-z population generation
from snapshot .2 and filtered halo indices, not established multi-epoch
stitching. Reconstructing that index from the public template still fails;
the exact local production mapping remains open. The matched-target model
comparison and pretraining experiment have not been run.

## Original next-test sequence

1. **Completed first controls.** Same current CG15 source, matched plot recipe,
   hard-label/probability comparison, then a common stellar-mass distribution.
   These substantially narrow the gap but leave different galaxy populations.
2. **Test selection within each sample.** In narrow redshift slices, derive a
   completeness limit using the actual targeting magnitude limits and mass
   convention; compare all available galaxies against the complete subset.
   Balance stellar mass and sky coverage. Do not extrapolate the 2025 low-z
   completeness spline to z=0.2–0.55 or call logM>9 complete there. Initially
   hold predictions fixed; rebuilding the tracer sample is a separate test.
3. **Check certainty within fixed samples.** The identical-row hard-label versus
   probability comparison above is complete. Next examine a high-
   confidence subset as a selection sensitivity, not a replacement population.
   More decisive old predictions can produce stronger-looking separation
   without being more accurate; probability weights can dilute associations.
4. **Check the physical scale and density dependence.** Recover the old label
   file's smoothing and threshold provenance before treating both sets of
   classes as the same environment. An experimental TNG script uses Rsmooth=2
   and threshold=0.2, whereas P12-A specifies 7 Mpc/h and 0.2; the experimental
   script alone does not prove the old checkpoint's label scale. Compare
   property trends at several directly measured neighbour scales and at fixed
   density. Current recorded P12-A/independent-density trends are already close.
5. **Then make a genuine model comparison.** Apply both methods to a shared,
   supported observation set or to an exposed mock with matched tracer
   selection, properties and environment definition. Simply running frozen
   P12-A below z=0.06 is outside its present range. Check performance against
   known mock environments, including any direct mass/colour inputs and the
   old feature normalization. Stronger property trends must not be the model-
   selection target, and a transductive simulation split does not establish
   performance on unseen survey volumes.

## Execution and verification

The audit, approved source scan and analysis ran with one CPU thread. No
inference, training, Slurm submission or release modification occurred.
Slurm allocation inspection failed on login13 with "Could not establish a
configuration source", including outside the sandbox. The user's explicit
bounded-login exception authorized this catalogue join and comparison.

Source: `workflows/presentation/legacy_gat_comparison_20261014.py`.
Use cosmic_env, unset PYTHONPATH/PYTHONHOME/PYTHONUSERBASE/LD_PRELOAD, and set
OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1, PYTHONNOUSERSITE=1.
`--stage audit --out OUTPUT` reproduces the completed work;
`--stage cigale --out OUTPUT` performs the join and seven-property
comparison, using the current supported P12-A z=.2–.3 wedge cache. The two
samples are balanced separately and cannot share a redshift distribution.
`--stage compare` reuses the joined old catalogue.

Plot exports were visually checked; ID/probability/selection checks passed.
All joined old coordinates/probabilities exactly reproduce the old pickle;
all retained IDs are unique integers. Fresh reads of 128 source rows reproduce
their saved properties exactly and pass position/redshift checks. An independent
Pandas grouped calculation reproduces every common-mass mean to <4e-15 and
both retained row counts exactly. The old pickle hash remains unchanged.
Full P12-A and CIGALE hashes are inherited from existing receipts, not freshly
recomputed; the old pickle and newly extracted outputs are freshly hashed.
`graphify` is unavailable; targeted inspection was used. No commit or push.
