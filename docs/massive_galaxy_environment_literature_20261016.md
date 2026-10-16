# Expected environmental trends above log10(Mstar/Msun)=10.69

Literature checked 2026-10-16. The threshold is 4.90e10 solar masses, not a
narrow mass bin. Our selected red population has mean log mass near 10.96;
the blue population is nearer 10.81. The appropriate expectation is weaker,
selection-dependent residual trends, especially within the blue population,
not a requirement that every property be flat or that red and blue have the
same slope. Published studies do not supply a transferable numerical prediction
for our CG15, observed-g−r-selected, probability-weighted P12-A statistic.

## What the literature actually measures

| Primary study | Relevant result | Limit on comparison with our analysis |
|---|---|---|
| [Peng et al. 2010, ApJ 721, 193](https://arxiv.org/pdf/1003.4747), Sections 3–4 | SDSS/zCOSMOS mass and density effects on population fractions are approximately separable; the star-forming population's sSFR is much less environment-dependent than its abundance. | Density and red fractions, not our inferred web classes or quenched fraction conditional on an observed-red cut. |
| [Wetzel, Tinker & Conroy 2012, MNRAS 424, 232](https://academic.oup.com/mnras/article/424/1/232/1008045), Figures 1–3 | At fixed stellar mass, active satellites and centrals have similar high-sSFR peaks, although satellites are more often quenched. The active peak does not show a significant halo-mass shift in their sample. | SDSS groups, median z≈.045, spectral SFRs and halo membership; not a bound on all cosmic-web effects or on CG15 SFR differences. |
| [Alpaslan et al. 2015, MNRAS 451, 3249](https://academic.oup.com/mnras/article/451/3/3249/1195624) | GAMA mass-matched populations have similar colours, energy outputs and morphologies across large-scale environments; unmatched samples exhibit stronger trends. | Mass conventions include h factors; environment and matching definitions differ. This is evidence against expecting large universal residuals, not an exact null prediction. |
| [Eardley et al. 2015, MNRAS 448, 3665](https://academic.oup.com/mnras/article/448/4/3665/976018) | GAMA tidal-web luminosity functions vary strongly, but their variation can be accounted for by local density in this analysis. | Luminosity functions, not SFRs. It demonstrates why a web association need not identify a geometrical effect beyond density. |
| [Kraljic et al. 2018, MNRAS 474, 547](https://academic.oup.com/mnras/article/474/1/547/4430643), Sections 4–5 | GAMA galaxies at .03≤z≤.25 show mass/type segregation around DisPerSE filaments and walls. At fixed mass, star-forming galaxies redden and lower sSFR near filaments. | A useful counterexample to a universal blue null: filament distance, rest-frame dust-corrected u−r selection and node exclusions differ from our four-class observed-colour comparison. |
| [Xu et al. 2020, MNRAS 498, 1839](https://arxiv.org/html/2009.07394), Sections 4.1–4.4 | In EAGLE centrals, after removing halo-mass effects, colour/sSFR web dependence weakens above 1.8e10 h^-1 Msun. Crucially, their Section 4.4 finds stronger high-mass trends when halo effects remain included. | Our mass/z controls do not remove halo mass or central/satellite mix. The halo-controlled null is therefore not the same estimand. This is simulation evidence at z=0, not a DESI prediction. |
| [Parente et al. 2024, ApJ 966, 154](https://arxiv.org/html/2312.10146v2), Sections V–VI | Their semi-analytic model produces much weaker sSFR/dust dependence above roughly 10^10.8 Msun; SDSS/GAMA comparisons broadly support convergence at high mass. The model threshold depends on its black-hole growth prescription. | Especially relevant mass scale, but not a universal boundary. Dust mass is not our fitted A_V. Their environments, sample, redshift and population statistics differ. |
| [Nandi & Pandey 2026, ApJ 1000, 251](https://arxiv.org/html/2507.18614v2), Sections II.5–III | SDSS DR18 shows sheet/cluster differences extending through log mass 10.5–11.5 using an 8-Mpc-smoothed tidal classification. | Their global density control restricts common density range rather than matching its distribution. Small-bin and coarse-bin sSFR contrasts differ substantially; no direct numerical calibration of our result follows. |

The newer SDSS study is a reason not to call high-mass environmental invariance
settled. Its selected high-mass/low-density bin has a quoted cluster median
log sSFR of −22.47 versus −10.80 in sheets, while the coarser version gives
−11.69 versus −10.99. I would audit SFR zeros/censoring and bin support before
using that extreme amplitude as a physical benchmark. That is a methodological
caution about transferability, not evidence that the paper's measurement is
incorrect. [Nandi & Pandey, Section III.2](https://arxiv.org/html/2507.18614v2).

## What we should expect for the current Loa statistic

1. **Quenched fractions and within-population SFRs are different quantities.**
   A substantial change in quenched fraction can coexist with a nearly fixed
   star-forming sequence. Our red cut still contains objects on both sides of
   sSFR=1e-11 yr^-1. Its mean log sSFR therefore mixes population abundance and
   changes within each star-formation state.
2. **A mass cut does not isolate a pure cosmic-web effect.** We control mass and
   redshift, but retain halo-mass, satellite-fraction and local-density differences
   between inferred environments. Those may contribute to the measured trend.
   A cluster-class probability is also not a measured galaxy-cluster membership.
3. **Blue-only changes may be weak, but not necessarily zero.** Conditioning on
   g−r<.75 selects surviving blue systems and can remove galaxies that redden as
   they quench. The observed cut is not a rest-frame star-forming definition;
   its population changes with redshift. More low-z blue galaxies need not mean
   an intrinsically stronger environmental effect.
4. **The remaining fitted-property signals need separate interpretation.** Age,
   A_V and FAGN are correlated outputs of SED fitting. FAGN is a fitted fractional
   contribution, not an AGN incidence rate; neither the sign nor size of a small
   fitted FAGN trend has a clean benchmark in the studies above.

A useful algebraic illustration, not a fitted explanation, is

    mean(log sSFR) = fQ * mean_Q(log sSFR) + (1-fQ) * mean_SF(log sSFR).

If the two state means were separated by 2 dex and did not change with environment,
a change in fQ of .086 would itself move the overall mean by −.172 dex.
Thus our measured −.208 dex within the colour-red sample cannot be interpreted
as suppression of the active sequence without decomposing the mixture.
This illustrative 2-dex separation was not measured or fitted in this analysis.

## Reading our measured amplitudes

At .20≤z<.22 and log mass≥10.69, the red controlled cluster-minus-void contrast
is +8.64 percentage points in fQ and −.208 dex in mean log sSFR. These establish
conditional associations with inferred environment under the recorded controls;
they do not establish a direct tidal/geometry mechanism. The blue contrast is
−.015 dex, conditional 95% [−.082,+.046], with only 1,462 retained objects.
A small blue signal is compatible with the literature and this uncertainty.
It is not evidence that model attenuation or domain transfer has been calibrated.

I would not assign a literature-derived expectation such as “exactly 5–10 points”
or “0.1 dex” to our statistic. A defensible amplitude test requires matching the
stellar-mass distribution, colour/SFR definition, halo or density controls,
redshift, environment scale/classification and probability-weighting procedure.
The most informative next interpretations are the lower-z rerun already requested,
a future matched-mass comparison between windows, and a decomposition by actual
star-formation state. The latter two are proposals, not analyses completed here.

Sources were checked as primary papers. No paper figures were republished, and
no local catalogues or code were uploaded to an external service for the review.
