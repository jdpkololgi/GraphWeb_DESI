# Cluster-deficit fix plan — graph-scale / N(z) (drafted 2026-06-19)

## Context (what we know)
The wedge NPE (and the Jraph regression) under-predict **clusters** on real DESI
(fraction ~0.027 vs Abacus truth ~0.058 at λ_th=0.2). Systematic falsification ruled
out under-density, Fingers-of-God, training-extrapolation and dense-structure shape
(see `SCIENCE_LOG.md` 2026-06-19). The one surviving, confirmed lever is a
**graph-scale domain shift**: DESI's wedge graph is 12% denser → edges 7% shorter →
scaled log(edge_length) sits **−0.11σ** below the Abacus-fit N(0,1) the GNN trained on
(`edge_scale.png`). It is shared by regression+NPE, so it lives in the GNN encoder.

**Key constraint discovered:** the n(z) mismatch is a *shape* difference — the mock is
slightly denser at z<0.22 but **sparser at z>0.23** (DESI/mock up to 1.36 near z≈0.27).
So the mock-side `--equal_data_dens y` / `--downsampling y` knobs (which *downsample*
the mock) go the **wrong direction** — they cannot raise the mock's high-z density to
match DESI. The mock fundamentally **under-produces high-z galaxies** (likely a
magnitude/k-correction or luminosity-function calibration in the CutSky BGS assignment,
since the r<19.5 cut then removes too many at high z).

## Phase 0 — cheap confirmation tests (NO retrain, NO mock rebuild)  [~1–3 h]
Goal: prove (or kill) graph-scale as the *cause* before spending days. Two tests,
cheapest first; success metric = **DESI cluster fraction moves toward 0.058** and the
edge-scale offset shrinks.

**0a. Inference-time edge-scale correction (~30 min, 1×A100).**
Add a `--edge-scale-shift` option to `infer_desi_wedge_flowjax.py`: after
`prepare_edges_for_jraph_forward`, add +0.111 to the scaled edge_length column (col 0)
so DESI's mean matches training (0.0). Re-run inference with the EXISTING linear model.
- If clusters recover → the −0.11σ edge shift is the cause and the fix is a trivial
  inference-time domain correction (no mock rebuild at all). Best case.
- Rebuild the diagnostics (`plot_lambda_th_sweep`, `eigenvalue_corner`,
  `plot_shape_misclassification`) on the corrected output.

**0b. DESI density-matching test (~2 h, if 0a is partial).**
Downsample the real DESI wedge to the mock's per-z-bin density, rebuild the DESI graph
(`build_desi_bgs_gudhi_graph` → `desi_graph_features_cugraph` → `subset_desi_graph_wedge`)
and re-run inference. This matches the graph SCALE end-to-end (not just the scalar
edge mean). If clusters recover here but not in 0a, the scale effect is structural
(connectivity), not just the edge-length mean → informs the Phase-1 route.

**Decision gate:** if neither 0a nor 0b moves the cluster fraction, graph-scale is NOT
the cause and we stop here (report the deficit as an open systematic; use λ_th=0.1).
If they do, proceed to a durable fix below.

## Phase 1 — durable fix (pick a route based on Phase 0)
### Route A (recommended, model-side, ~1 day): scale-invariant edge features + retrain
Make the GNN robust to overall graph density so it transfers regardless of n(z).
1. In `build_abacus_sbi_cache.py` edge construction (and the DESI parity helper
   `abacus_gnn_parity.py`), replace absolute `log(edge_length)` with a **per-graph
   normalized** edge length (e.g. `log(edge_length / median_edge_length)` or a
   per-graph z-score), so the feature is invariant to global density. Keep direction
   and density_contrast as-is. Add a `--scale-invariant-edges` flag so we can A/B it.
2. Rebuild the Abacus linear SBI cache (`--linear-increments --power-scale-node-features
   --scale-invariant-edges`). Validate (round-trip, splits) as before.
3. Retrain the linear NPE: `jraph_sbi_flowjax.py --increment_mode linear --epochs 7000`,
   4×A100 (~3 h). Re-run the eval (TARP/SBC/class-fractions) — confirm the Abacus-side
   result is unchanged (still well-calibrated, R²≈0.82, cluster ≈0.056).
4. Re-run DESI inference with the matching scale-invariant edge transform. The DESI
   wedge graph is unchanged.

### Route B (sim-realism, slower, only if a faithful mock is wanted): densify the mock
The mock under-produces high-z galaxies, so the fix is to ADD galaxies, not downsample.
This is non-trivial and `--equal_data_dens` will NOT do it. Steps:
1. Diagnose the high-z deficit in `upstream_prepare_mocks_Y3_bright.py` /
   the CutSky magnitude assignment: compare mock vs DESI N(z) and the r-mag
   distribution by z; check whether the r<19.5 cut + mock k-corrections remove too
   many high-z galaxies (likely culprit).
2. Re-tune the magnitude/luminosity assignment (or HOD) so the mock passes the right
   number at z>0.23, then re-run **prepare → fiberassign → mkCat → maglim**
   (fiberassign is the multi-hour step; reuse it only if geometry is unchanged).
3. Rebuild graph/cache, retrain, re-infer (as Route A steps 2–4).
This is the most faithful fix but the riskiest to land before the talk — recommend it
as post-conference work unless Phase 0 shows scale is the *whole* story.

## Phase 2 — graph + cache + train + infer (for the chosen route)  [~half day compute]
1. **Mock catalog** (Route B only): new `mock_bgs_maglim.fits`.
2. **Graph-ready FITS:** `write_fiberassign_mock_science_fits.py`.
3. **Graph:** `build_abacus_graph.py --no-apply-y1y5-filter --no-exclude-invalid-box-index`.
4. **Features (rapids-gnn, GPU):** `abacus_graph_features_cugraph.py`.
5. **Wedge subset:** `subset_abacus_graph_wedge_for_sbi.py` + `subset_cugraph_metrics_for_wedge.py`.
6. **SBI cache:** `build_abacus_sbi_cache.py --three-targets-only --linear-increments
   --power-scale-node-features [--scale-invariant-edges]`.
7. **Train:** `jraph_sbi_flowjax.py --increment_mode linear --epochs 7000` (4×A100).
8. **Eval:** `plot_flowjax_posteriors.py` (TARP/class-fractions) — Abacus side unchanged.
9. **DESI inference:** `infer_desi_wedge_flowjax.py` with the new model/cache/edge npz;
   then the full diagnostic suite (`plot_lambda_th_sweep`, `eigenvalue_corner`,
   `plot_shape_misclassification`, `plot_edge_scale`).

## Verification / success metrics
- **Primary:** DESI inferred cluster fraction recovers toward Abacus truth (0.058),
  and the dense inferred-cluster rate at matched density/shape (the 92% residual in
  `plot_shape_misclassification.py`) shrinks.
- **Edge-scale:** `plot_edge_scale.py` DESI scaled-edge mean → ~0 (was −0.11σ).
- **No regression on Abacus:** retrained model keeps TARP calibration, R²≈0.82,
  Abacus cluster ≈0.056, ordering-violation ≈0.
- **Sanity:** class fractions for void/wall/filament stay matched (don't break the
  parts that already transfer well).

## Recommended sequencing before the talk
1. Phase 0a (30 min) — likely decisive and possibly a complete fix.
2. If partial, Phase 0b (2 h).
3. If scale confirmed and time allows: Route A retrain (~1 day) for a durable result.
4. Route B (mock densification) = post-conference.
Fallback for the deck regardless: present the falsification chain + edge-scale finding
as a characterized systematic, with λ_th=0.1 as the pragmatic recovery knob.

## Key paths
- Mock scripts: `/pscratch/sd/d/dkololgi/abacus/SecondGen_Mocks/ph000/scripts/`
  (`upstream_prepare_mocks_Y3_bright.py --downsampling`,
  `upstream_mkCat_SecondGen_amtl.py --equal_data_dens`).
- Abacus wedge npz (edge scaler): `.../graph_constructions/wedges/path1_fiberassign/..._cugraph_gnn_arrays.npz`.
- Linear cache: `.../abacus/sbi_caches/path1_flowjax_3d_lineareig/...`.
- Linear model: `.../abacus/sbi_runs/path1_wedge_flowjax_3d_Bcorrected_linear/`.
- DESI wedge: `.../graphweb_desi/outputs/desi_wedge_expanded_..._from_fullgraph/`.
- Inference + diagnostics: `GraphWeb_DESI/workflows/sbi_inference/`.
