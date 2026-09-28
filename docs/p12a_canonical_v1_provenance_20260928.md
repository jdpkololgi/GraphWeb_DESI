# Canonical v1 generator investigation — 2026-09-28

## New source lead

Canonical means our comparison path, not an independently certified production recipe.
The sole file in `/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v1/z0.200/`
is owned by jpiat and dated 2024-11-27 (7175868480 bytes). It is therefore not an
independent producer route around Jade, contrary to the implication of the recent
recommendation. We can nevertheless investigate publicly available code without a reply.

Located https://github.com/JadePiat/hodpy with four public branches. Master is
identical to amjsmith/hodpy master; inspecting master alone misses the relevant code.
The abacus branch is pinned at 14222dcf794b3bff68987bf8e271203db20f8730.
Its head is dated 2025-02-04, AFTER the canonical file timestamp. The latest
ancestor returned before 2024-11-28 is a2473ad00a46ab75cbbaef646671997727f8682d
(2024-11-21). Neither is an established executed version.

## Concrete differences from the old shared_code generator

Source inspection of the pinned abacus branch:

- `make_mocks/make_catalogue_cutsky_abacus.py` provides the snapshot-to-cutsky
  entrypoint, defaults PHOTSYS=S, snapshot z=0.2, snapshot magnitude floor -18,
  faint apparent limit 20.2, and separate low-z unresolved/resolved populations.
- `hodpy/lookup.py` selects `hod_fits_c000_ph000_wsys.txt` and a luminosity
  function predicted by the fitted HOD. Its comment explicitly aims to preserve
  clustering at z=0.2. Both wsys/nowsys tables are present; their names alone do
  not establish which observations/weights were used in fitting.
- DESI-specific N/S K-correction, colour distribution and observed/rest-colour
  tables replace the older GAMA/SDSS mapping. These tables are included in the
  public tree, not merely named in a missing private configuration.
- The E-correction default is Q=0.67, pivot z=0.1. The BGS LF is described as
  already E-corrected; the HOD LF uses P=Q=0 where instantiated. This is NOT the
  same prescription as the older evolving LF P=1.8,Q=0.7. These parameters
  occur at different stages and cannot be compared as a one-number change.
- `cut_sky_tools.py` assigns DESI colours separately to centrals/satellites,
  computes apparent magnitudes through DESI_KCorrection, then applies the flux
  cut. This is a specific, testable explanation for failure of the old colour
  identity on v1, not proof yet that it reproduces that file.
- Satellite radial extent defaults to 1.5 R200; this and the HOD changes mean
  matching n(z) does not guarantee matching galaxy-matter relationships.

## Additional branches and cautions

`my_abacus` (0b555378a2a854cdf57310f950259608ff0f6f88) changes luminosity distance
from (1+z_obs)*chi(z_obs) to (1+z_obs)^2/(1+z_cosmo)*chi(z_cosmo) by default,
and changes central velocities. The shown central-velocity patch assigns an
all-galaxy-length draw to a satellite-only slice: it needs validation before use.
`fix_velocity` (5894fbbafa7f4b2e4d69ad37852d38d528cb2aef) instead correctly sizes
that draw to Nsat. These are candidate explanations of prior velocity differences,
not identification of which branch generated either catalogue.

## Next executable tests

1. Compare the November-2024 ancestor to the pinned branch, including lookup
   table hashes; establish historically possible candidates.
2. On already exposed ph000, replay N/S DESI colour and apparent-magnitude
   mappings against canonical v1. Test both distance conventions separately.
   Use deterministic rows and residual distributions, not just total counts.
3. If a mapping matches, compare HOD/LF tables and retained magnitude-redshift
   distributions; reconstruct one controlled phase using that candidate recipe.
4. Apply the documented Loa selection/assignment chain and evaluate clustering
   and galaxy-matter closure before any retraining or VAC replacement.

No new numerical catalogue test, generation, correction or retraining occurred
in this source-provenance pass. The exact executed recipe remains unresolved.
GitHub API branch/tree/diff responses are retained under
`evidence/p12a_canonical_v1_provenance_20260928/`.

## Alex upstream comparison and numerical replay (2026-09-28)

Alex's https://github.com/amjsmith/hodpy/tree/abacus resolves to EXACTLY the
same commit as Jade's abacus branch: 14222dcf794b3bff68987bf8e271203db20f8730.
The newer BGS prescription is upstream Alex code, not a Jade-specific change.
GitHub comparison of pre-file ancestor a2473ad to that head changes only
`lookup/uchuu/central_fraction_uchuu.npy`; the BGS code/tables are unchanged.
Branch and comparison responses are preserved in
`evidence/p12a_canonical_v1_replay_20260928/`.

Ran `workflows/p12a_vac/replay_v1_photometry.py` on CPU59007827/nid004176.
Every1000th row of the 74,748,522-row ph000 file, then .15<=Z<.55, gives58,825
rows spanning the file (not a sky-prefix sample). No new phases exposed.
Bundled BSD-licensed upstream K-correction source and N/S tables are hashed in
RESULTS.json. Colour tests are independent of cosmology; magnitude distances
use approximate flat LCDM c000 in Mpc/h, not an exact CLASS replay.

Single-region mappings fail globally: colour RMS .04122N/.02467S mag.
Selecting the lower colour residual per row gives13,283N/45,542S, with
colour RMS9.61e-8mag, max3.04e-6mag. Those labels separate at DEC32.375:
N minimum32.37818, S maximum32.37301. This motivated an explicit fixed-declination
replay check, recorded separately in RESULTS.json rather than treating per-row
fitted choices as an independently specified selection.

With these labels, apparent magnitude RMS5.95e-5mag, max1.03e-4mag for the
upstream observed-redshift luminosity distance and Q=.67 prescription.
The24,195 sampled r<19.5 rows have colour RMS5.14e-8mag and magnitude
RMS5.10e-5mag. Small magnitude residuals are consistent with the approximate
distance calculation, but that explanation is not yet a verified exact replay.

Conclusion: strong numerical identification of canonical-v1 DESI N/S
photometric mapping in Alex's abacus branch. This does NOT establish the executed
HOD fit, random seed, complete generation history, Loa selection parity, or
which change accounts quantitatively for the high-z population gain.
Next: isolate HOD/LF versus photometric effects with fixed-parent comparisons;
verify exact distance and original N/S assembly before regenerating a suite.

Final fixed DEC>32.375 check has ZERO label disagreements and identical residuals.
The alternative mixed z_cosmo/z_obs luminosity distance yields .01022mag RMS
(max .09278), substantially worse than upstream z_obs. All three bounded
steps completed0 (20s/19s/19s); allocation released. Finite-value and shell-count
assertions passed. No fitting, regeneration or VAC changes.
