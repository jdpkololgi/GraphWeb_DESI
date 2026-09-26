"""Independent serialized-product QA for the provisional Loa canary."""
import argparse, hashlib, json
from pathlib import Path
import fitsio
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): h.update(block)
    return h.hexdigest()

def main(root):
    complete=json.loads((root/'VAC_COMPLETE.json').read_text())
    assert sha(complete['catalogue'])==complete['catalogue_sha256']
    vac=fitsio.read(complete['catalogue']); ids=vac['TARGETID']
    assert len(ids)==complete['rows'] and np.all(np.diff(ids)>0)
    assert np.all((vac['Z']>=.15)&(vac['Z']<.55))
    assert np.all((vac['QUALITY']&32)!=0)
    covered=[]
    for item in complete['draw_shards']:
        assert sha(item['path'])==item['sha256']
        with np.load(item['path']) as shard:
            target=shard['TARGETID']; draws=shard['eigenvalue_draws']; ix=np.searchsorted(ids,target); rows=vac[ix]
            assert np.array_equal(ids[ix],target) and draws.shape==(len(target),512,3)
            supported=rows['SUPPORTED']; d=draws[supported]
            assert np.isfinite(d).all() and np.all(np.diff(d,axis=-1)>=0)
            assert np.isnan(draws[~supported]).all()
            assert np.array_equal((rows['QUALITY']&1)==0,supported)
            if supported.any():
                assert np.array_equal(rows['EIGENVALUE_MEAN'][supported],d.mean(axis=1))
                q=np.quantile(d,[.05,.16,.5,.84,.95],axis=1).astype('f4')
                for j,name in enumerate(['Q05','Q16','Q50','Q84','Q95']):
                    assert np.array_equal(rows['EIGENVALUE_'+name][supported],q[j])
                classes=(d>.2).sum(axis=-1)
                prob=np.stack([(classes==i).mean(axis=1) for i in range(4)],axis=1).astype('f4')
                assert np.array_equal(rows['P_WEB'][supported],prob)
                assert np.allclose(prob.sum(axis=1),1)
            for name in ['EIGENVALUE_MEAN','EIGENVALUE_Q05','EIGENVALUE_Q16','EIGENVALUE_Q50','EIGENVALUE_Q84','EIGENVALUE_Q95','P_WEB']:
                assert np.isnan(rows[name][~supported]).all()
            covered.extend(target.tolist())
    assert np.array_equal(np.sort(covered),ids)
    provenance=json.loads((root/'SOURCE_MANIFEST.json').read_text())
    for name,value in provenance['observation_source_sha256'].items():
        assert sha(root/'source_snapshot'/name)==value
    candidate=Path(provenance['illustris_source']).parent
    assert sha(candidate/'RUN_MANIFEST.json')==provenance['illustris_manifest_sha256']
    manifest=json.loads((candidate/'RUN_MANIFEST.json').read_text())
    for name,value in manifest['source_sha256'].items():
        assert sha(candidate/'source'/name)==value
    supported=vac['SUPPORTED']
    result=dict(status='PASS',rows=len(vac),supported_rows=int(supported.sum()),unsupported_rows=int((~supported).sum()),unique_sorted_targetids=True,joint_draws_per_supported_row=512,ordered_finite_draws=True,saved_summaries_and_class_probabilities_exact=True,all_provisional=True,source_files_verified=len(manifest['source_sha256'])+len(provenance['observation_source_sha256']),quality_bit_counts={str(bit):int(((vac['QUALITY']&bit)!=0).sum()) for bit in [1,2,4,8,16,32]},mean_class_probabilities=vac['P_WEB'][supported].mean(axis=0).tolist(),science_release_ready=False,completion_sha256=sha(root/'VAC_COMPLETE.json'))
    out=root/'VAC_VALIDATION.json'
    with out.open('x') as f: json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);main(p.parse_args().root)
