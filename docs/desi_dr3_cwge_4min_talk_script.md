# DESI DR3 KP — Cosmic Web & Galaxy Environment parallel session
## 4-min talk — slide script v2 (2026-07-07; updated for session schedule)

**Fixed talk title (from schedule):** *Learning the cosmic web: dynamical
environment value added catalogues for DESI BGS* — slide 1 must use this.

**Session context (order matters — you are talk 5 of 7):**
1. Tojeiro — DR2 DisPerSE catalogues + mock validation
2. Li — MTV (mass/tidal/velocity) field reconstruction from BGS groups ← nearest
   methodological neighbour (also targets the tidal field)
3. Torres-Gomez — ASTRA cosmic-web classification VAC ← another classification VAC
4. Rincon — void environment & low-z galaxies (DR2 BGS) ← science consumer of
   environment labels, immediately before you
5. **You**
6. Lee — DisPerSE persistent-homology/topology VACs
7. Kim — cosmology + neutrino mass with filaments

Consequences:
- **Zero cosmic-web setup.** By talk 5 the room has seen T-web-style
  classifications and DisPerSE several times. Spend the saved time on the
  differentiators.
- **Positioning, not introduction:** the deck's job is "how this VAC differs from
  the ones you just saw" — (a) *dynamical* target (tidal tensor eigenvalues, the
  quantity in the title), (b) *probabilistic*: full per-galaxy posterior, hence
  threshold-free class probabilities and propagatable uncertainties, (c)
  *calibration-tested* on mocks before touching data (TARP/SBC), (d) truth-free
  closure on the data itself.
- **Bridge lines to neighbours** (cheap goodwill, saves words): slide 2 —
  "unlike the geometric and classification catalogues you've just seen, the
  target here is the tidal tensor itself"; slide 5 — "for exactly the science
  Hernan just showed, these labels come with error bars you can propagate."
- **No per-talk Q&A** — consolidated Q&A (10 min, all speakers) + 40-min
  discussion. So: do NOT bank on Q&A to complete an argument; every claim must
  land inside the 4 minutes. The Q&A prep below doubles as discussion prep.

**Audience:** DESI-internal, environment-fluent. No cosmic-web motivation needed,
no NPE/flow tutorial needed. They care about: what the environment measure IS,
whether it's trustworthy on real BGS data, and what environmental science it enables.

**Throughline:** *"A calibrated posterior over the tidal tensor for every galaxy
in the trained DESI BGS wedge — and the known environmental trends fall straight
out of it."* (SCOPE: wedge subvolume only — trained on the matched Abacus
subvolume/z-shell; the full BGS footprint is out-of-distribution and would need
retraining. Do NOT claim a whole-catalogue VAC.)

**Budget:** 5 slides, ~45–55 s each. Numbers below are reconciled to the final
figures (SI production run `desi_wedge_flowjax_linear_si`, closure join 2026-06-21,
spotlight figures 2026-07-06).

Figure root (NERSC): `/pscratch/sd/d/dkololgi/graphweb_desi/figures/`

---

### Slide 1 — Title (0:00–0:15)
**Title (fixed by schedule):** Learning the cosmic web: dynamical environment
value added catalogues for DESI BGS
**Subtitle:** per-galaxy tidal-tensor posteriors — graph neural network + neural
posterior estimation, trained on AbacusSummit, applied to Loa
**Visual:** `desi2026_spotlight/sky_desi_inferred.png` (full-bleed background) +
`desi2026_spotlight/legend_classes.png`
**Say (2 sentences):** "This is a DESI BGS wedge where every galaxy carries an
*inferred* cosmic-web environment — not a hard label, a posterior. I'll show you
how we get it, why you can trust it, and what it's already doing for environmental
science."

---

### Slide 2 — Dynamical definition → why posteriors (0:15–1:00)
**Title:** Dynamical definition of environment (keeps "dynamical" from the title)
**Hero visual = the collapse LADDER** (the draft-slide-28 ladder, kept — the
arrows are physically correct: N eigenvalues above λ_th = collapse along N axes,
Void 0 → Wall 1 → Filament 2 → Cluster 3; this is the Zel'dovich/tidal-web
picture and is exactly why the target is "dynamical"). Skewer video REMOVED; the
pivot line below does its discrete→continuous job.
**Layout (equation + ladder + 3 short lines only — cut the 4 λ-inequality rows,
the 4 side notes, and the eigenvalue histogram → backup):**
- top corner, small: `T_ij = ∂²Φ/∂x_i∂x_j → λ₁ ≥ λ₂ ≥ λ₃`
- anchor line over the ladder (caption for the arrows): "Each eigenvalue above
  λ_th = one axis of gravitational collapse"
- caption under the ladder: "λ_th = 'collapsed within a Hubble time' — 0
  (Hahn+07) vs 0.44 (Forero-Romero+09): arbitrary"
- pivot line at the bottom (replaces skewer video + 'discrete imperfect' notes):
  "Real galaxies sit between the rungs → we infer the full posterior over
  (λ₁,λ₂,λ₃), not a hard label"
**Opening bridge line:** "You've just seen geometric and classification
catalogues — here the target is the tidal tensor itself, per galaxy, as a
posterior." (Tidal tensor = the same object Qingyang Li reconstructs — free nod.)
**Say:** the ladder is the discrete picture; the pivot line is the whole talk in
one sentence — thresholding posterior samples then gives calibrated class probs
at *any* λ_th, with error bars. One clause on magnitude mattering physically
(tidal-torque/IA couples to the tensor, not the class).
**Physics caveat (for you, not the slide):** eigenvalue>0 = "collapses eventually"
(linear theory); λ_th>0 = collapse completes within a Hubble time (the FR09
justification). Keep arrows INSIDE each rung (collapse directions of that
environment) — do NOT let a between-row arrow imply voids evolve into clusters;
voids are a distinct (expanding) fate, not an early cluster stage.

---

### Slide 3 — How: mock-trained, recovers the true distributions (1:00–1:50)
**Title:** Trained on AbacusSummit BGS — recovers the true eigenvalue distributions
**Pipeline strip on the RIGHT (text/diagram, reuse Cambridge pipeline slide):**
AbacusSummit CutSky BGS mock (fiberassign + n(z)-matched to Loa) → Delaunay graph
on observed redshift-space positions → GNN encoder → normalizing-flow posterior
(NPE) per galaxy
**Bullets:**
- Trained purely on the mock (~100k-galaxy fiberassign wedge); applied to DESI
  with per-graph feature normalisation to absorb the mock-vs-data graph-scale shift
- **Recovered eigenvalue distributions track the truth** (grey = Abacus truth,
  blue = Abacus NPE, magenta = DESI NPE overlaid): the trained-vs-true agreement
  is the "it works" evidence for THIS audience — skip TARP/SBC (right idea, wrong
  room; they read distribution alignment, not coverage curves)
- Zero-shot on DESI: class fractions land on the Abacus truth — V/W/F/C
  0.23/0.43/0.30/0.05 vs truth 0.27/0.41/0.26/0.06 — including the rare 5–6%
  cluster class
**Visual/layout:** pipeline diagram RIGHT, VERTICAL 3-row eigenvalue distribution
LEFT — `desi_wedge_flowjax_linear_si/eig_dist_vertical_3way.pdf` (all three
series; use `eig_dist_vertical_abacus_only.pdf` to hold DESI back until slide 4).
Vector PDF, native 7.0×7.3 in, labels/ticks 25 pt — place at native size for
25 pt on the slide; script `workflows/sbi_inference/plot_eig_dist_vertical.py`
(edit FS to rescale). `class_fractions_comparison.png` can move to slide 4 or sit
here small.
**Say:** "the recovered eigenvalue distributions land on the truth — and transfer
to DESI without retraining." One breath on per-graph normalisation ("the mock and
the data disagree about absolute graph scale; the model sees only relative
geometry"). TARP/SBC calibration stays in backup for anyone who asks.

---

### Slide 4 — On DESI: cosmic-web cartography + closure (1:50–3:00) — HERO
**Title:** Applied to DESI Loa: the trends you know fall out of the inferred web
**Visual:** `desi_wedge_cigale_hz/spotlight_cartography_closure.png` (CIGALE-HZ
version — regen'd 2026-07-09) — single composite: Abacus training wedge (true
T-web) | DESI BGS wedge (NPE-inferred) sky maps on top; quenched fraction / (g−r)
/ sSFR vs E[Σλ] closure panels below. This one figure IS the slide.
**Say:**
- "~111k BGS Bright galaxies in this Loa wedge, each with a posterior; ~100k with
  CIGALE-HZ SED masses+SFRs (90% match) for the closure."
- "Truth-free closure test: against the *inferred* environment, quenched fraction
  rises 0.39→0.58 void→cluster, galaxies redden 0.78→0.91 in (g−r), median sSFR
  falls −10.2→−11.6 (1.4 dex) — monotonic in class AND continuous in the
  tidal-tensor trace (ρ_s = +0.11 / +0.14 / −0.14), at N≈100k significance."
- Honest clause: "environment is inferred from positions, so a density–quenching
  correlation is partly expected — the point is the *T-web decomposition* carries
  the signal, with calibrated uncertainty."

---

### Slide 5 — Mass–colour by environment + take-home (3:00–4:00)
**Title:** Blue cloud → red sequence, void → cluster — and you can use this
**Visual (preferred):** `desi_wedge_cigale_hz/sfms_by_environment.png` — with the
CIGALE-HZ SFRs the SFR–M* plane is now CLEANLY BIMODAL (blue cloud + red sequence
+ green-valley gap), and the blob balance visibly shifts star-forming→quenched
Void→Cluster (the FastSpecFit version smeared this out). `mstar_gr_by_environment.png`
(2×2 M*–(g−r) contours) is the alt: Void bimodal, blue lobe shrinks to a single
red-sequence peak by Cluster. Median (g−r) rises 0.780→0.825→0.868→0.907; median
log M* (CIGALE) 10.72→10.77→10.81→10.85.
**Alt visuals:** `desi_wedge_cigale_hz/{sfms_all,mstar_gr_all_bimodal}.png` (whole-
wedge classic bimodality), `mstar_gr_overlay.png` (class-coloured contours),
`sfms_eigen_continuous.png` (SFR–M* coloured by mean λ). Scripts
`plot_sfms_environment.py`, `plot_mstar_color_environment.py` (both take
`WEDGE_PARQUET`/`WEDGE_FIGDIR`); properties = CIGALE-HZ re-join `desi_wedge_cigale_hz`.
**Data note:** full wedge N=111,503 (110,251 unique); ~100,679 (90%) have CIGALE-HZ
SFR/mass. The ~22k figure = the thin 3° Dec fan slice, NOT the whole wedge.
**Bullets:**
- The colour bimodality shifts blue→red with environment — the red/quenched
  population grows void→wall→filament→cluster (consistent with the closure trend,
  shown here as the full mass–colour distribution not just a fraction)
- Deliverable: per-galaxy (λ₁,λ₂,λ₃) posteriors + class probabilities **within
  the trained wedge** (RA 120–160°, z 0.20–0.30). These relations are validated
  for wedge galaxies only — the model is trained on the matched Abacus subvolume,
  so the full BGS footprint (lower z, different n(z)) is OUT of distribution and
  would need RETRAINING on footprint-spanning mocks. Not a whole-catalogue VAC yet.
- Use cases for this room: environmental quenching at fixed mass, IA/spin–tidal
  alignment (the tensor, not the class), void catalogues with membership
  probabilities, threshold-free comparisons between web-finders
**Bridge line (Rincon spoke immediately before):** "for exactly the void–galaxy
science we just saw, these labels come with uncertainties you can propagate."
**Say (closing line):** "If your science needs 'which environment is this galaxy
in?' — this gives you that answer *with error bars*, validated on mocks and
closing on the data."

---

## Consolidated Q&A + discussion prep

There is NO per-talk Q&A: 10 min consolidated Q&A across all 7 speakers (expect
1–2 questions max aimed at you), then 40 min discussion. Two consequences:

- **Consolidated Q&A:** questions will often be comparative ("how do these VACs
  relate?"). The one to expect: **"How does this compare to Qingyang Li's MTV
  tidal-field reconstruction / to ASTRA / to DisPerSE?"** Pocket answer:
  complementary object — reconstructions and web-finders give a *deterministic*
  field or label from a chosen pipeline (group finder, persistence threshold,
  λ_th); this VAC gives the *posterior* over the tensor per galaxy, so it can (a)
  reproduce any of their labelings threshold-free and (b) put error bars under
  exactly the science Hernan showed. Cross-comparison on a common wedge is the
  natural KP exercise.
- **Discussion (40 min) — points worth seeding:**
  1. Propose a **common validation wedge**: all environment VACs (DisPerSE ×2,
     ASTRA, MTV, this) evaluated on one footprint with a shared mock — the
     Abacus fiberassign BGS wedge machinery already exists and is offered.
  2. Offer the **closure test as a shared, truth-free validation metric**: any
     VAC can be scored on whether f_quench/(g−r)/sSFR trends emerge against its
     environments at fixed M* — pipeline runs in ~1 min per catalogue given
     TARGETIDs (`build_desi_wedge_property_join.py` is `--preds`-flagged).
  3. **Uncertainty propagation as a KP standard**: downstream environmental
     statistics should carry environment-label uncertainty; posterior samples
     make that mechanical.

### Likely questions (answers)

1. **Fibre incompleteness?** Mildly environment-dependent in Loa; the mock is
   fiberassign-processed and n(z)-harmonised to match, so the model sees the same
   selection. Biases absolute fractions at most, not the trends; the closure trend
   is robust.
2. **RSD / Fingers-of-God?** Everything (training + inference) is in observed
   redshift space — the model learns environment *in the presence of* RSD rather
   than pretending it away. FoG was explicitly ruled out as the driver of the
   early cluster deficit; that was a graph-scale (n(z)) domain shift, fixed by
   per-graph normalisation.
3. **Why a GNN and not DTFE/classical density + tidal solve?** Measured classical
   floor on the identical test split: DTFE + exact FFT tidal solve reaches
   R² = 0.55/0.64/0.66 (λ₁/λ₂/λ₃) vs 0.78/0.81/0.89 for the GNN+NPE — the learned
   headroom is real, largest for λ₁ and λ₃.
4. **Smoothing scale of the target field?** ~7 Mpc/h Gaussian; tested 5–15 —
   global optimum is broad, clusters favour the fine end; 7 is the adopted
   compromise.
5. **vs DisPerSE / other web-finders?** Different object: they give geometry/hard
   labels; we give the posterior over the underlying tensor — you can *derive*
   consistent labels at any threshold and propagate uncertainty into downstream
   statistics. A quantitative cross-comparison is natural DR3 KP work.
6. **Graph construction sensitivity?** Delaunay on both training and inference
   sides (no mismatch); alpha-complex tried and worse. Edge features are
   pairwise geometry; per-graph normalisation removes absolute-scale dependence.
7. **When/where can I get the labels? / Can you run the whole BGS catalogue?**
   Wedge catalogue (parquet: TARGETID, posterior summaries, class probabilities)
   exists NOW for the trained subvolume (RA 120–160°, z 0.20–0.30). It is NOT a
   whole-footprint VAC: the model is trained on the matched Abacus subvolume, so
   the rest of BGS (lower z, different n(z)/density regime) is out-of-distribution
   — a full-catalogue VAC needs RETRAINING on footprint-spanning mocks (amortised
   inference makes *running* cheap, but does not license extrapolation off the
   training volume). Honest framing: a validated wedge demonstrator + a clear path
   to the VAC, not the VAC itself yet.
8. **Is the closure circular?** Partly by construction for *density*; the
   non-trivial content is (a) the T-web class decomposition ordering, (b)
   survival at fixed M*, (c) calibrated widths — a mislabeled or miscalibrated
   web would wash these out (and did, before the SI fix, in the cluster class).

## Backup slides (have, don't present)
- `desi_wedge_flowjax_linear_si/tarp_linear.png` + `sbc_linear.png` (calibration)
- `desi_wedge_flowjax_linear_si/eig_dist_buildup_with_desi.png` (domain-shift story)
- `desi_wedge_flowjax_linear_si/cluster_recovery_bars.png` + `edge_scale.png`
  (cluster-deficit diagnosis + fix)
- `desi_wedge_flowjax_linear_si_closure/closure_continuous_{trace,p_cluster}.png`
- `desi_wedge_flowjax_linear_si/posterior_width_sky_map.png`, `class_sky_map.png`
- `desi2026_spotlight/closure_{quenched,gr,ssfr}.png` (single-panel closure, if a
  question wants one metric big)

## Deck review v3 (2026-07-07, vs ASTRA deck; draft slides 26–32)

**What ASTRA's deck does right (emulate):** one idea per slide; one big visual
per slide with almost no competing text; a "wow" full-bleed fan plot with zero
words; coloured keyword highlighting instead of dense bullets; citations as tiny
footnotes. Their content is *thinner* than ours — the polish is layout economy.

**Slide-by-slide on the current draft:**
- **27 (wedge video):** drop the video. Playback risk + it eats ~20 s of
  attention for one message. Replace with the static **fan pair** (new
  `fan/fan_abacus_truth` + `fan/fan_desi_inferred`, generation in progress):
  trained-on (truth) → applied-to (inferred), full-bleed, legend only. This
  becomes the deck's wow slide AND absorbs slide 30's sky-map panel.
- **28 (dynamical definition):** REVISED — KEEP the collapse ladder (it's
  physically correct and is the "dynamical" hook; see Slide 2 above for the
  decluttered layout). The clutter to cut is the FOUR λ-inequality text rows, the
  FOUR side notes, and the eigenvalue histogram — NOT the ladder/arrows. Reduce to
  equation + ladder + anchor line + λ_th caption + pivot line. Skewer video
  removed; the pivot line carries discrete→continuous.
- **29 (How / train on Abacus):** keep the pipeline graphic; delete the three
  eigenvalue histograms (calibration detail → backup). One proof line: "TARP
  coverage at the MC noise floor before touching data". Check wording "HOD
  mocks" — the production training set is the CutSky fiberassign+n(z)-matched
  wedge; say what's true.
- **30 (Inference on DESI):** currently ~4 visuals (sky map, repeated GraphNET
  diagram, fractions inset, histograms). One message: zero-shot transfer,
  fractions land on truth. Hero = `class_fractions_comparison.png` alone; the
  sky map moves to the fan slide; never repeat the pipeline diagram.
- **31 (closure):** right content. Use `closure_strip.png` (or the categorical
  version) + the LEFT panel of `closure_mass_control.png` — "and it's not just
  mass" is the beat this audience scores you on.
- **32 (summary):** fine at 3 bullets; availability line must be SCOPED — "wedge
  demonstrator now; full-footprint VAC needs retraining on footprint-spanning
  mocks," NOT "whole-catalogue VAC." Don't imply the model runs on all of BGS.

**Net structure (5 content slides):** title → fan pair (trained→applied) →
lean definition/why-posteriors → method+calibration+fractions → closure+mass →
summary. Matches the v2 timing plan above; the fan pair replaces both videos.

**Fan-plot technique (for reference):** matplotlib polar axes — theta = RA
(deg→rad) with `set_thetamin/thetamax` clipped to the wedge so it renders as a
partial fan, r = z (or comoving distance, tick-labelled in z), thin ~2–4° Dec
slice so structure stays filamentary, tiny points, class colours on black.
DESI's famous slice visuals (Claire Lamman) are exactly this construction;
`desiutil.plots` covers all-sky projections but not fans — plain matplotlib is
the tool. Script + settings will live in `figures/desi2026_spotlight/fan/`.

## Timing notes
- Hard cut discipline: slides 4–5 are the payload for this room; if running long,
  compress slide 2 (the "why posteriors" argument can survive as one spoken
  sentence over slide 1's map).
- The composite hero (slide 4) needs ~15 s of silent look-time — build the
  sentence rhythm around letting them read the three closure panels.
