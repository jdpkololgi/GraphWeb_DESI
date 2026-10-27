"""Representative CPU audit on a copy of preflight data; originals unchanged."""
import os
from pathlib import Path
import tarfile
import time
import h5py
import numpy as np
import loa_survey_production_20261016 as p
from audit_loa_survey_preview_20261016 import audit_tile, normalizer_scale
r=p.r


def main():
    r.require_compute(); p.check_campaign(p.ROOT)
    root=p.ROOT/'RECOVERY_VALIDATION'; tile='NGC_005_008_003'
    # This new archive is made only from the already validated recovery copy.
    p.packing.pack(root,tile,4)
    row=next(x for x in r.read(root/'INDEX.json')['cases'] if x['id']==tile)
    # NumPy warns that HDF5 dtype annotations are not serialized in NPZ.
    # Verify actual catalogue values independently on a fresh prepared tile.
    annotations={}
    with h5py.File(p.ROOT/'catalogue.h5','r') as catalogue, np.load(p.ROOT/'tiles'/tile/'prepared/arrays.npz') as saved:
        selection=slice(row['row_start'],row['row_start']+row['row_count'])
        for name in catalogue:
            original=catalogue[name][selection]
            np.testing.assert_array_equal(saved[name],original)
            if original.dtype.metadata:
                annotations[name]=str(original.dtype.metadata)
    start=time.monotonic(); report,ids=audit_tile(root,row,normalizer_scale(root))
    for model in report['models'].values():
        if sum(model['retained_class_draw_counts'])!=4*len(ids):
            raise ValueError('Retained class counts do not reconcile')
        if sum(model['eligible_class_draw_counts'])!=4*report['eligible']:
            raise ValueError('Eligible class counts do not reconcile')
    names=('audit_loa_survey_preview_20261016.py','test_loa_preview_audit_20261016.py',
           'validate_loa_preview_audit_20261016.py','loa_survey_preview_audit_20261016.slurm')
    sources={str(Path(__file__).with_name(name)):r.sha(Path(__file__).with_name(name)) for name in names}
    destination=p.ROOT/'PREVIEW_AUDITOR_SOURCES.tar'
    with tarfile.open(destination,'x') as archive:
        for path in sources: archive.add(path,arcname=Path(path).name,recursive=False)
    result=dict(passed=True,tile=tile,fine_draws=24,retained=len(ids),eligible=report['eligible'],
        seconds=time.monotonic()-start,sources=sources,source_archive_sha256=r.sha(destination),
        exact_archive_checks=True,array_physics_checks=True,class_count_checks=True,
        stored_coarse_draws_checked=8,saved_core_mass_reconciled=True,
        fresh_prepared_catalogue_values_exact=True,dtype_annotations_not_serialized=annotations,
        scope='one dense recovered preflight tile, not a complete survey preview audit',
        original_preflight_preserved=True,job=os.environ['SLURM_JOB_ID'])
    r.atomic_json(p.ROOT/'PREVIEW_AUDITOR_VALIDATED.json',result)
    print(result,flush=True)


if __name__=='__main__': main()
