# Loa observed-colour comparison at 0.15 ≤ z < 0.20

Completed 2026-10-16 using approved CPU allocation 59857020 on nid200019; allocation released and final squeue query returned no job or stderr.

## Selection and support

The colour-audited fixed limit is log10(Mstar/Msun) ≥ 10.60 (3.9811e10 Msun), only 0.09 dex below the earlier 0.20–0.22 limit of 10.69. Of 534,070 selected galaxies, 518,962 are red, 14,909 blue and 199 outside the colour cuts. The blue count is 7.76 times the earlier 1,922.

The same Milky-Way-corrected observed-colour cuts apply: red 0.75 < g−r < 1.5; blue g−r < 0.75. There is no K correction. Hard void/sheet/filament/cluster counts are red 98,173/222,675/167,396/30,718 and blue 3,607/6,733/4,036/533. Every class exceeds 100, but the predeclared 20,000-per-colour gate still fails. This remains exploratory; the count rule is not a power calculation.

## Controlled results

Controls use a common mass/redshift distribution across environments separately within each colour, with 0.1 dex mass and 0.01 redshift cells; the mass panel controls redshift only. Soft controls retain 518,374 red and 14,123 blue galaxies for the non-mass properties. Below are cluster-minus-void contrasts with conditional 95% sky-bootstrap intervals (256 draws). Quenched contrasts are percentage points; logarithmic quantities are dex.

| Property | Red contrast [95% CI] | Blue contrast [95% CI] |
|---|---:|---:|
| quenched | 9.37050 [8.93756, 9.74141] | 0.61767 [0.02508, 1.21450] |
| logmass | 0.04089 [0.03917, 0.04305] | -0.00108 [-0.00725, 0.00352] |
| logssfr | -0.20174 [-0.21207, -0.19058] | -0.02307 [-0.04245, -0.00206] |
| logsfr | -0.20139 [-0.21183, -0.19037] | -0.02394 [-0.04323, -0.00282] |
| colour | 0.02750 [0.02529, 0.02972] | 0.00425 [0.00160, 0.00659] |
| logage | 0.01115 [0.01034, 0.01188] | 0.00136 [-0.00088, 0.00381] |
| av | -0.05622 [-0.05990, -0.05167] | -0.00426 [-0.01694, 0.00751] |
| fagn | 0.00355 [0.00275, 0.00427] | 0.00010 [-0.00300, 0.00341] |

Red quenched fractions rise from 62.90% in voids to 72.27% in clusters after these controls. The red log-sSFR contrast is −0.202 dex; these associations persist with finer and joint soft/hard controls. A red-sample mean includes changes in the quenched/star-forming mixture and is not a fitted star-forming sequence.

The blue log-sSFR contrast is small: −0.023 dex (95% CI −0.042 to −0.002), with a +0.00425 mag colour contrast. Blue quenched fraction differs by +0.618 percentage points, but its interval includes zero with finer and joint controls. At 0.05 dex with joint soft/hard support only 6,451 blue objects remain, and both log-sSFR and colour intervals include zero. These are weak, control-sensitive hints across multiple properties, not a robust detection across all specifications. No blue age, attenuation or AGN-fraction contrast is established.

## Interpretation and validation

Changing redshift also changes surveyed volume and observed-colour membership: this is not an isolated causal test of stellar mass. The intervals condition on fixed selection and environment predictions; they exclude SED errors, completeness uncertainty, selection systematics and simulation-to-data transfer. They do not establish calibration or causation. See [massive-galaxy literature review](massive_galaxy_environment_literature_20261016.md) for physical context.

Independent checks passed for all ten colour percentiles, all 128 raw/controlled mean vectors, 256 valid bootstrap draws in each calculation, FITS roundtrip and fresh original-VAC parity for 257 sampled rows across all original columns. The fixed threshold exceeds every old bin limit, proving the cached source selection does not truncate qualifying rows. Both final PNGs were visually inspected. Small archived outputs have verified SHA256 parity; source/output hashes are in the manifests.

## Products

- [Property curves](figures/loa_colour_lowz_20261016/properties_z015020.png) and [PDF](figures/loa_colour_lowz_20261016/properties_z015020.pdf).
- [Contrasts](figures/loa_colour_lowz_20261016/contrasts_z015020.png) and [PDF](figures/loa_colour_lowz_20261016/contrasts_z015020.pdf).
- [Evidence and CSV](evidence/loa_colour_lowz_20261016/ARCHIVE_INDEX.json).
- Full preserved-column FITS: `/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_colour_lowz_20261016_v1/LOA_COLOUR_MASS_COMPLETE_z015020.fits`; large data remain external to Git.
