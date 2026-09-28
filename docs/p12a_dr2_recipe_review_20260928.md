# DR2 wiki, CAI slides and mock pipeline review — 2026-09-28

## Assessment

These sources document related BGS problems and materially strengthen the
observational-selection investigation. They do not establish that our old
Abacus raw-parent high-z deficit is already fixed in a named release.
Separate intrinsic-parent differences from assignment completeness: missing
competitors can change observed n(z), but cannot explain a difference already
present in a raw parent before assignment.

## Slide-deck findings

Reviewed text of the 309-page CAI Updates2026 deck, exported PDF, and visually
checked relevant BGS slide36 plus image-only slides29-30 (LRG reddening
conclusion/blank). Source:
https://docs.google.com/presentation/d/1xLu19oZcjZELENla4LMvagxZKjQS_uQ7T2GV5qqtxUw/edit

- Slide36, September3: BGS redshift-dependent completeness mismatch, suspected
  incomplete full BGS population and missing low-z cluster targets. Figure is
  HOLI/GLAM BGS-21.35 at approximately0.1<z<0.4, not our flux-limited sample to
  z=.55. This is a related documented mechanism, not identification of our bug.
- Slides128-130, May21: complete/on-the-fly BRIGHT samples compared against
  BGS_ANY-02 altMTL with a redshift-dependent absolute cut did not match; reported
  density difference~5%. Matching the BRIGHT definitions was proposed. This
  affects AbacusSummitBGS_v2/kibo-v1 products directly, but concerns a restricted
  absolute-magnitude sample; it is not proof of a50% full-BRIGHT correction.
- Slide215, April9: low-z clipping and possible compensating density boost in
  HOLI/GLAM proposed for assignment statistics. Later slide36 warns that the
  resulting completeness can still differ. Counts alone are not validation.
- BGS status text discusses testing high-priority white dwarfs as
  competitors for a~3% completeness excess. Exact executed recipe not established.
- September24 slides4-9 describe DR3/Nevis preparation, time-dependent bad
  fibres, DR11 target retirement and incorporating main targets observed on1b
  tiles. Matterhorn matching requires release-specific history, not a name swap.
- Slides13-29 are image-only; slide29 is an LRG-reddening conclusion that
  explicitly proposes BGS as a future extension. These pages were not all
  individually inspected; the BGS-specific text and slide36 were reviewed.

## Attached wiki: what the versions mean

The supplied page is the DR2 **data LSS catalogue** release history. Loa v1.1,
v2 and v2.1 must not be conflated with raw Abacus BGS v1 or mock-pipeline versions.

- Loa v1: redshifts differ from Kibo at~0.1% level; one additional BRIGHT tile;
  page says no mock rerun needed for that transition. Kibo here is DA2/DR2,
  not evidence that we used a DR1 sky by virtue of its mountain name.
- Loa v1.1: FRAC_TLOBS_TILES bug fixed. Essential completeness audit, not an
  asserted fix to intrinsic mock galaxy abundance. Current comparison usesv2.1.
- Loa v2: time-dependent bad-fibre/petal masks and BGS clustering range expanded
  to0.002<z<0.6. Expanding output range does not create high-z mock galaxies.
- Loa v2.1: corrected altMTL replication affecting PIP completeness weights;
  page states other catalogue types unaffected. Our full_HPmapcut comparison
  is not interchangeable with a PIP clustering catalogue.
- BGS WEIGHT_ZFAIL performance explicitly flagged as inadequate in the wishlist.
  We must retain assignment-only versus assignment*zfail sensitivity and inspect
  the linked TSNR/redshift validation; current weights are not parent truth.
- Official *_full_HPmapcut_nz.txt files contain weighted Nbin and header area:
  useful independent normalization check with matched samples/units/areas.

## GitHub audit at LSS d942b990860e016e7558c740966fea3e5ce7e7f3

| Source | Finding | Disposition |
|---|---|---|
| Sandbox/mockAMTL_DA2_LSSpipe.txt | Explicit BRIGHT-specific preparation, imaging mask, ledgers, altMTL, potential assignments, catalogue stages | Retain staged recipe and separate success markers |
| prepare_mocks_Y3.py | Unconditional top-level exit after target-bit dictionary; incomplete BGS branches and inconsistent arguments/targets variable | Not a drop-in replacement for our pinned BRIGHT adapter |
| prepare_mocks_Y3_bright.py | Existing source of our adapter; includes faint competitor prescription | Audit full competitor support; don't assume generic DARK examples supersede it |
| test_target_density.py | Projected-density/area comparison helper; private paths and potential shared random-output writes | Adapt read-only matched-area diagnostic, do not execute defaults |
| prepare_holi_bgs.py + calibrate_nz_prep.py | Explicit Loa ANY and BRIGHT n(z) matching; WISE fraction, faint HIP, tail construction; some R_MAG_ABS=-21/-19 placeholders | Confirms real n(z) inputs exist; assignment/covariance recipe, not faithful photometric training labels |
| ab2ndgen_bgsbright02loa_interactive.sh | BRIGHT full then BRIGHT-02 redshift-dependent absolute cut, usepota, Loa | Preserve BRIGHT full sample; don't inherit -02 science selection |
| glam_bgsv2_DA3DA2_lss_pipeline_interactive.sh | survey=DA3 but surveycat=DA2 and specdata=loa-v1 | DA3 provenance does not imply Matterhorn output or DR3 footprint |
| linked mock_Y3 wrappers | DARK/EZmock examples; some use CPU shared allocation and stage-specific memory; LSS --nz produces n(z)/FKP bookkeeping | Operational examples, not BGS generator calibrations or universal resource settings |
| getpotaDA2_mock.py/readwrite_pixel_bitmask_da2.py | Potential assignments and additional mask stages | Potential-assignment contract needed before LSS completeness; LRG-specific mask not a BGS cut |
| fiberassign PR471 | Reuse original fibreassign stuck-sky status to avoid repeated sky lookup | Performance change, not luminosity/population repair |
| brickmask desi branch | MPI brick-mask assignment | Performance option; require rowwise mask/NOBS parity before replacing current implementation |

Important distinction: upstream n(z)-calibration randomly retains rows according
to target-density / parent-density. It cannot supply missing galaxies where
that ratio exceeds1. A surplus permits such selection but a1D match does not
establish colour, luminosity, clustering or galaxy-matter consistency. It is
legitimate to reproduce this documented recipe as a labelled control, distinct
from determining the physical/selection origin of the discrepancy.

Source references:
https://github.com/desihub/LSS/tree/d942b990860e016e7558c740966fea3e5ce7e7f3/scripts/mock_tools
https://github.com/desihub/LSS/blob/d942b990860e016e7558c740966fea3e5ce7e7f3/Sandbox/mockAMTL_DA2_LSSpipe.txt
https://github.com/desihub/fiberassign/pull/471
https://github.com/cheng-zhao/brickmask/tree/desi
Source snapshots, static assertions and hashes: evidence/p12a_dr2_recipe_review_20260928.

## Revised priority, with no scientific repair applied

1. Audit BRIGHT/FAINT/HIP/WISE and other competing targets across the entire
   target redshift range, including below our science z=.15. Compare angular
   density versus local density, tile coverage and magnitude. Our adapter includes
   BRIGHT+FAINT but uses random69.5% FAINT retention/20% HIP and MWS_TARGET=0;
   it has no explicit WISE branch. Account for any real standards/secondary inputs
   before concluding that all stellar competitors are absent.
2. Enforce the same BRIGHT parent and magnitude definition for complete,
   assigned, vetoed and selected comparisons. Test the documented ANY/-02 mismatch
   against exact files, not just filenames or shared cuts.
3. Check official weighted Nbin/areas and redshift-success diagnostics with
   matched cuts, retaining our assignment-only bracket. v2.1 fixes do not make
   the weighted catalogue an imaging-complete population.
4. Complete ph000 forward replay as an explicitly approximate competitor model;
   measure residual completeness by z/magnitude/local density. Existing queued
   preparation59013807 remains unchanged, not silently retuned by this review.
5. Keep the recovered photometry/LF tests: raw-parent differences cannot be
   explained entirely by downstream assignment. Independent LF conventions and
   clustering/galaxy-matter gates remain necessary.

No catalogues altered, no additional jobs submitted by this review. Batch59013807
was still pending when checked. No proof of a ready Loa- or Matterhorn-matching
replacement, and no VAC science-release promotion.
