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
