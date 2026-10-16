"""Independent frozen-input/derived-product validation after rendering."""
import os
import socket
import json
from pathlib import Path
import numpy as np
import fitsio
from run import OUT, BASE, WINDOWS, _previous, masks, sha, write


def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Approved compute allocation required')
    windows=json.loads((OUT/'windows.json').read_text())
    origins={'z015020':BASE/'loa_colour_lowz_20261016_v1',
             'z020022':BASE/'loa_colour_split_20261016_v2',
             'z020030':BASE/'loa_colour_split_20261016_v1'}
    checks=[]
    for key,lo,hi,oldlim,expected in WINDOWS:
        old=fitsio.read(origins[key]/f'LOA_COLOUR_MASS_COMPLETE_{key}.fits',columns=['TARGETID'])
        now=fitsio.read(OUT/f'observed_original_{key}.fits',columns=['TARGETID'])
        assert np.array_equal(np.sort(old['TARGETID']),np.sort(now['TARGETID']))
        for variant in ['rest_primary','paired']:
            path=OUT/f'{variant}_{key}.fits';a=fitsio.read(path)
            assert np.unique(a['TARGETID']).size==len(a)
            assert ((a['Z']>=lo)&(a['Z']<hi)&a['REST_VALID']).all()
            limit=windows[key]['rest_limit' if variant=='rest_primary' else 'paired_limit']
            assert (np.log10(a['MASS_CG_15'])>=limit).all()
            assert np.array_equal(a['GR_REST01'],a['ABSMAG01_SDSS_G']-a['ABSMAG01_SDSS_R'])
            ix=np.linspace(0,len(a)-1,min(257,len(a)),dtype=int)
            orig=fitsio.read(_previous.VAC,rows=a['SOURCE_ROW'][ix])
            for col in orig.dtype.names:
                assert np.array_equal(orig[col],a[col][ix],equal_nan=True) if orig[col].dtype.kind in 'fc' else np.array_equal(orig[col],a[col][ix]),col
            checks.append(dict(window=key,variant=variant,rows=len(a),fresh_original_vac_rows=len(ix),source_columns=len(orig.dtype.names)))
    vectors=0;raw_boot=[];control_boot=[];unavailable=0
    for p in OUT.glob('relationships_*.json'):
        results=json.loads(p.read_text())
        for props in results.values():
            for rec in props.values():
                for v in rec['raw'].values():
                    if v is None:continue
                    raw_boot.append(v['bootstrap_valid']);vectors+=1
                    assert np.isclose(v['mean'][-1]-v['mean'][0],v['contrast'],atol=1e-12)
                for control in rec['controlled'].values():
                    if control['status']!='measured':unavailable+=1;continue
                    control_boot.append(control['bootstrap_valid'])
                    assert np.isclose(sum(control['reference_weights']),1,atol=1e-12)
                    for v in control['methods'].values():
                        vectors+=1;assert np.isclose(v['mean'][-1]-v['mean'][0],v['contrast'],atol=1e-12)
                        assert np.isfinite(np.array(v['mean'])).all()
    assert min(raw_boot+control_boot)>=205
    for manifest,folder in [('prepare_manifest.json',OUT),('analysis_manifest.json',OUT),('presentation_manifest.json',OUT/'presentation')]:
        d=json.loads((OUT/manifest).read_text())
        for name,h in d['outputs'].items():assert sha(folder/name)==h,(manifest,name)
    result=dict(job=os.environ['SLURM_JOB_ID'],host=socket.gethostname(),code_sha256=sha(Path(__file__)),
        source_checks=checks,curve_vectors_checked=vectors,unavailable_controls=unavailable,
        raw_bootstrap_valid_range=[min(raw_boot),max(raw_boot)],controlled_bootstrap_valid_range=[min(control_boot),max(control_boot)],
        original_observed_targetid_sets='Exact agreement for all three windows',all_manifest_hashes='passed')
    write('validation.json',result);print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
