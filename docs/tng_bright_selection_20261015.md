# TNG z=0 approximate BGS Bright selection at Loa observing redshifts

The selection-only experiment substantially reduces the quenched-fraction
offset, but its residual is sensitive to the photometric approximation.
No new snapshots, model inference, training, graph construction, or Slurm jobs.

| Shared population | TNG minus Loa before brightness cut | After approximate colour-corrected cut | Reduction |
|---|---:|---:|---:|
| All supported mass–redshift cells | 29.32 pp | 4.70 pp | 24.61 pp |
| 10.5 ≤ log10(M/Msun) < 11.5 | 28.99 pp | 7.23 pp | 21.76 pp |

Offsets are means across four equal-weight environment quarters. Every arm
uses identical observed Loa mass–redshift cell weights. Conditional spatial
95% intervals for the residual are [3.99,5.44] and [6.31,8.16] pp;
paired reduction intervals are [24.03,25.22] and [21.13,22.34] pp.
These intervals exclude photometric, SFR, mass, and transfer systematics.

![Selection overlay](/global/u2/d/dkololgi/fmsu_tokyo_20261010/tng_bright_selection_20261015/selection_overlay.png)

## Inputs and construction

- Existing held-out TNG300-1 z=0 catalogue and R7 truth eigenvalues from
  `~/fmsu_tokyo_20261010/tng_density_control_20261015`; original unique IDs,
  stellar masses, instantaneous total gas SFRs, and spatial blocks retained.
- Read g/r `SubhaloStellarPhotometrics` from 600 existing snapshot-99 group
  files, joining by global Subfind ID. All selected cached g-r values agree.
  These are summed stellar photometry, not simulated survey measurements.
  [Official TNG specifications](https://www.tng-project.org/data/docs/specifications/)
  document the field and band order.
- Same geometrically supported Loa wedge and CG15 eligibility as the previous
  overlay: exact TARGETID alignment to targeting cache, 0.20≤z<0.30,
  mass 9≤logM<11.5, explicit r<19.5. All wedge targeting limits are 19.5;
  no northern-limit approximation is needed here.
- Place the same TNG population at ten observing-redshift bin midpoints
  0.205,...,0.295. Use Planck18 luminosity distances; do not evolve the
  galaxies, duplicate them into an independent volume, or alter environments.
- Approximate observed r as M_r + distance modulus + K(z), with
  K(z)=(g-r) ln(1+z)/ln(6231/4770) - 2.5 log10(1+z).
  This is a local power-law f_nu approximation between rest-frame SDSS g/r,
  evaluated at the redshifted r-band wavelength. It does not integrate
  spectra through Legacy Survey filters or model dust, emission lines,
  spectral breaks, apertures, surface brightness, noise, or fiber selection.
- Sensitivities: K=0, and colour K with uniform +0.3-mag dimming. The latter
  is an arbitrary sensitivity test, **not** a dust model or fitted correction.
- Require ≥80 galaxies in every arm per 0.25-dex mass × 0.01-redshift cell.
  Common support comprises 58 cells/87,256 Loa galaxies, 97.0% of 89,948
  eligible rows. High-mass support includes all 40 cells/68,918 eligible
  galaxies. Low-mass support is redshift dependent: do not read the left
  panel as a complete rectangular mass range.
- Reweight each arm to those same Loa cell counts; matching is binned, not
  exact continuous mass/redshift matching. No SFR or quenching weights.
- Rank TNG R7 trace and Loa posterior-mean R7 trace separately within cells,
  recomputing quarters after selection; ties have fractional membership.
  Thus quarters are equal ranks, not proven identical physical environments.
  Quenched means sSFR<1e-11/yr, with original instantaneous versus CG15
  measurement mismatch preserved.
- 200 spatial bootstrap draws: TNG 5×5×5 cube blocks, Loa original sky
  regions; all repeated observing-redshift placements share the same TNG
  draws. Cell weights and support are fixed. No independent-copy inflation.

## Interpretation and sensitivities

Colour-corrected selected TNG fractions run 51.1–56.1% across quarters,
versus Loa 44.9–51.8%; high-mass TNG 64.4–70.5% versus Loa 55.3–63.6%.
Environmental contrasts change from 12.67 to 4.93 pp, versus Loa 6.89 pp;
high-mass 7.88 to 6.10 pp, versus Loa 8.30 pp. Selection improves absolute
agreement but does not establish agreement of the environmental slope.

At the same reference cell weights, the approximate cut retains 40.0% of
quenched versus 81.7% of star-forming TNG galaxies; high-mass retentions
are 50.2% versus 92.6%. This directly demonstrates selection on luminosity
at fixed coarse mass bins, beyond the earlier mass-mixture correction.

The mean residual is +15.55 pp with no K correction and -6.28 pp with
additional 0.3-mag dimming (high mass +20.14 and -6.38 pp). These are not
a calibrated confidence envelope. They demonstrate that the apparent
near-agreement depends strongly on luminosity modelling. Do not report
the nominal ~84% gap reduction as an identified fraction of the real
observational discrepancy explained by selection.

Keeping the preselection parent rank boundaries gives selected TNG
50.4,52.6,52.8,56.6% (high mass 63.5,66.1,66.5,71.1%), close to reranked
curves. The large vertical shift is therefore not solely rank redefinition.

Next discriminating inputs are observer-frame, dust-attenuated synthetic
photometry and comparable SFR measurement, followed by true epoch evolution.
This experiment does not validate P12-A amplitudes or measure real-DESI truth.

## Reproducibility and validation

Code: `workflows/presentation/tng_bright_selection_20261015.py`.
Run with cosmic_env Python; bounded single-thread catalogue/plot diagnostic
under the existing user-authorized login fallback.
Outputs: `~/fmsu_tokyo_20261010/tng_bright_selection_20261015/` contains PNG/PDF,
`curves.csv`, `RESULTS.json`, source-bound `photometry.npz`, original chunk
metadata, and `validation.json`.

348 direct cell-mean checks passed. Eight independent scalar quenched-mean
aggregations reproduced the averages of the rank curves using the equivalent
spectral-index K formula; high-mass uncut TNG and Loa curves reproduce the
previous overlay (10 checks total, max discrepancy 2.22e-16).
Colour parity, ID alignment, nested dimming selection, finite bootstrap
denominators, and weight normalization passed. PNG inspected visually.
Graphify is unavailable (command not found); targeted direct inspection used.
