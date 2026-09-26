# Does the Loa sample contain unwanted galaxies?

Reviewed existing complete-catalogue censuses and historical targeting replay;
no new catalogue scan, selection change or inference performed.

All5,436,413 VAC science TARGETIDs match original desitarget1.1.1 bright inputs
and pass the complete historical BRIGHT replay (TARGET_REPLAY.json). This verifies
target eligibility, not that every measured redshift or galaxy classification is
correct. Earlier nominal component failures are superseded by that full replay.
Every selected row passes BGS bit2, GOODHARDLOC, GOODPRI, LOCATION_ASSIGNED and
positive NOBS_G/R/Z in the full census. Existing input is full_HPmapcut, not raw
unvetoed targeting. Exact historical mask/version equivalence to mocks remains open.

A real quality-definition difference exists: current GraphWeb uses ZWARN=0,
DELTACHI2>=25, SPECTYPE=GALAXY. Adding DELTACHI2>40 while keeping GALAXY removes
13,293/5,436,413 rows (0.2445%), including1,853/137,328 (1.349%) at0.45–0.55.
These lower-confidence rows are not all demonstrated bad redshifts. The earlier
quality40 replacement comparison also dropped the GALAXY condition and could
increase counts; it must not be confused with this additional exclusion.

The full census finds80 selected rows with nonzero COADD_FIBERSTATUS. They need
bit-level interpretation against the Loa production's accepted-status rules;
nonzero alone has not been established as equivalent to an invalid observation.

Northern PHOTSYS targeting legitimately extends to nominal extinction-corrected
r<19.54, versus19.5 south. Uniform r<19.5 removes about5.3% of high-z NGC and
almost none in SGC, but changes the intended sample. It is an alignment control,
not proof of wrongly admitted targets. The full historical replay includes SGA
recovery and Gaia rejection; blanket component cuts can reject valid recoveries.
All baseline common-sky galaxies pass the tested TSNR2_BGS>1000 condition.

Recommendation: test DELTACHI2>40 retaining GALAXY as the conservative quality
variant, resolve the80 fibre-status rows, and pin the Loa mask/observation recipe.
Propagate any adopted sample change through counts, response fields and mock
observation modelling. Do not add an absolute-magnitude/BRIGHT-02 cut unless
that different scientific sample is intended. The tested stricter conjunction
still leaves a large high-z discrepancy; none establishes missing cuts as its cause.

Sources: docs/evidence/p12a_loa_full_20260925/SELECTION_AUDIT_V2.json;
docs/evidence/p12a_targeting_alignment_20260925.json;
docs/evidence/p12a_closure_20260926/TARGET_REPLAY.json.
