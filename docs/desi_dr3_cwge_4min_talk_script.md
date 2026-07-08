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

**Throughline:** *"A calibrated posterior over the tidal tensor for every BGS galaxy —
and the known environmental trends fall straight out of it."*

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

### Slide 2 — What we infer, and why posteriors (0:15–1:00)
**Title:** The T-web, but with uncertainties
**Opening bridge line:** "You've just seen geometric and classification
catalogues — here the target is the tidal tensor itself, per galaxy, as a
posterior."
**Bullets (max 3):**
- Target = the 3 eigenvalues (λ₁, λ₂, λ₃) of the tidal tensor at each galaxy —
  the *continuous* quantity behind every T-web classification
- Hard V/W/F/C labels need an arbitrary λ_th (0 in Hahn+07, 0.44 in
  Forero-Romero+09); a posterior over (λ₁,λ₂,λ₃) subsumes them all — thresholding
  posterior samples gives calibrated class probabilities at *any* λ_th
- Per-galaxy uncertainty: boundary galaxies get honestly broad posteriors instead
  of a confident wrong label
**Visual:** `desi_wedge_flowjax_linear_si/eigenvalue_corner_environments.png`
(corner coloured by class regions), or a single skewer frame from the Cambridge
deck if it reads better at a glance.
**Say:** stress "your favourite threshold, with error bars" — that's the hook for
this room. One clause on magnitude mattering physically (tidal-torque/IA couples
to the tensor, not the class).

---

### Slide 3 — How: mock-trained, calibration-tested (1:00–1:50)
**Title:** Trained on AbacusSummit BGS, calibrated before touching data
**Pipeline strip (text/diagram, reuse Cambridge pipeline slide):**
AbacusSummit CutSky BGS mock (fiberassign + n(z)-matched to Loa) → Delaunay graph
on observed redshift-space positions → GNN encoder → normalizing-flow posterior
(NPE) per galaxy
**Bullets:**
- Trained purely on the mock (~100k-galaxy fiberassign wedge); applied to DESI
  with per-graph feature normalisation to absorb the mock-vs-data graph-scale shift
- Calibration on held-out mock: TARP coverage within the Monte-Carlo noise floor
  (max |ECP−α| = 0.015, N=3000) + flat SBC ranks
- Zero-shot on DESI: class fractions land on the Abacus truth — V/W/F/C
  0.23/0.43/0.30/0.05 vs truth 0.27/0.41/0.26/0.06 — including the rare 5–6%
  cluster class
**Visual:** `desi_wedge_flowjax_linear_si/class_fractions_comparison.png`
(hero, right half) + small `tarp_linear.png` inset (or TARP to backup).
**Say:** "We do not show you anything on data that wasn't calibration-tested on
the mock first." One breath on per-graph normalisation ("the mock and the data
disagree about absolute graph scale; the model sees only relative geometry").

---

### Slide 4 — On DESI: cosmic-web cartography + closure (1:50–3:00) — HERO
**Title:** Applied to DESI Loa: the trends you know fall out of the inferred web
**Visual:** `desi2026_spotlight/spotlight_cartography_closure.png` — single
composite: Abacus training wedge (true T-web) → arrow → DESI BGS wedge
(NPE-inferred) on top; quenched fraction / (g−r) / sSFR vs E[Σλ] closure panels
below. This one figure IS the slide.
**Say:**
- "111,000 BGS Bright galaxies in this Loa wedge, each with a posterior. Joined
  to FastSpecFit at 99.7%."
- "Truth-free closure test: against the *inferred* environment, quenched fraction
  rises 0.74→0.84 void→cluster, galaxies redden 0.78→0.91 in (g−r), median sSFR
  drops half a dex — monotonic in class AND continuous in the tidal-tensor trace,
  at N=111k significance."
- Honest clause: "environment is inferred from positions, so a density–quenching
  correlation is partly expected — the point is the *T-web decomposition* carries
  the signal, with calibrated uncertainty."

---

### Slide 5 — Not just mass + take-home (3:00–4:00)
**Title:** It's not just stellar mass — and you can use this
**Visual:** `desi_wedge_flowjax_linear_si_closure/closure_mass_control.png`
(left panel is enough: f_quenched vs class in 4 log M* bins, all rising).
**Bullets:**
- The environment trend survives at fixed M*: quenched fraction rises
  void→cluster within *every* stellar-mass bin (e.g. 0.88→0.94 for
  log M* > 11) — not a re-derivation of the mass–environment relation
- Deliverable: per-galaxy (λ₁,λ₂,λ₃) posteriors + class probabilities for BGS —
  wedge today, footprint-scalable (amortised: inference is a forward pass) —
  planned as a VAC
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
7. **When/where can I get the labels?** Wedge catalogue (parquet: TARGETID,
   posterior summaries, class probabilities) exists now; footprint runs are
   amortised (no retraining) — VAC intended. (Adjust promise level to comfort.)
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

## Timing notes
- Hard cut discipline: slides 4–5 are the payload for this room; if running long,
  compress slide 2 (the "why posteriors" argument can survive as one spoken
  sentence over slide 1's map).
- The composite hero (slide 4) needs ~15 s of silent look-time — build the
  sentence rhythm around letting them read the three closure panels.
