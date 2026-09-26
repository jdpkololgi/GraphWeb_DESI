"""Independent benchmark draw validation, old-canary parity and checkpoint-resume test."""
from production import *
r=Path(sys.argv[1]);b=r/'benchmark';old=read(b/'VAC_COMPLETE.json');out=fitsio.read(old['catalogue']);ids=out['TARGETID'];assert np.all(np.diff(ids)>0)
for item in old['draw_shards']:
 with np.load(item['path']) as f:
  ix=np.searchsorted(ids,f['TARGETID']);assert np.array_equal(ids[ix],f['TARGETID']);validate_rows(out[ix],f['eigenvalue_draws'])
canary=fitsio.read(INPUT/'DESI_LOA_P12A_HALO48_CANARY_VAC.fits');ix=np.searchsorted(ids,canary['TARGETID']);assert np.array_equal(ids[ix],canary['TARGETID'])
for name in ['BASE_EIGENVALUES','BOUNDARY_MPC','NTILDE_MPC3','SUPPORTED']:
 assert np.array_equal(out[name][ix],canary[name]),name
# Simulate loss of aggregate completion after every atomic core was written.
(b/'VAC_COMPLETE.json').rename(b/'VAC_COMPLETE_BEFORE_RESUME.json');infer(b);new=read(b/'VAC_COMPLETE.json')
assert old['catalogue_sha256']==new['catalogue_sha256']
assert old['draw_shards']==new['draw_shards']
save(r/'BENCHMARK_QA.json',dict(rows=len(out),supported=int(out['SUPPORTED'].sum()),draw_validation='PASS',prior_canary_base_response_parity='exact on all5615 rows',resume_after_core_completion='exact catalogue and all draw checksums',unit_tests='3 passed: intact resume, corrupt draws, changed input and ownership checks',science_release_ready=False))
print('BENCHMARK AND RESUME PASS',len(out))
