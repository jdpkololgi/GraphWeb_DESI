"""Write portable curve tables and a final output manifest on the compute node."""
import csv
import json
import os
from pathlib import Path
import socket
from run import OUT, VAC, PHOT, SED, sha, write_json


def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Requires allocated compute node')
    expected={'validation.json':124,'relationship_validation.json':12,'relationship_validation_full.json':12}
    for name,count in expected.items():
        assert json.loads((OUT/name).read_text())['passed']==count
    fields=['scope','selection','mass_bin_dex','property','method','environment','mean','ci95_low','ci95_high','matched_rows','valid_rows']
    contrasts=['scope','selection','mass_bin_dex','property','method','cluster_minus_void','ci95_low','ci95_high','matched_rows']
    with (OUT/'curves.csv').open('x',newline='') as f,(OUT/'contrasts.csv').open('x',newline='') as g:
        w=csv.writer(f);w.writerow(fields);v=csv.writer(g);v.writerow(contrasts)
        for scope,name in [('primary_interval','relationships.json'),('full_selected_range','relationships_full.json')]:
            for selection,widths in json.loads((OUT/name).read_text()).items():
                for width,props in widths.items():
                    for prop,record in props.items():
                        if record['status']!='measured':continue
                        for method,values in record['methods'].items():
                            for j,env in enumerate(['void','wall','filament','cluster']):
                                w.writerow([scope,selection,width,prop,method,env,values['mean'][j],values['ci95'][0][j],values['ci95'][1][j],record['rows'],record['input_rows']])
                            v.writerow([scope,selection,width,prop,method,values['contrast'],*values['contrast_ci95'],record['rows']])
    with (OUT/'mass_limits.csv').open('x',newline='') as f:
        w=csv.writer(f);w.writerow(['zlo','zhi','population','rows','faint_rows','p90','p95','adopted_logmass95'])
        for b in json.loads((OUT/'mass_limits.json').read_text()):
            for population,r in b['populations'].items():
                w.writerow([b['zlo'],b['zhi'],population,r['rows'],r.get('faint_rows'),r.get('p90'),r.get('p95'),b['logmass_limit']])
    source=Path(__file__).parent
    manifest={'status':'complete','job_id':os.environ['SLURM_JOB_ID'],'node':socket.gethostname(),
              'sources':{str(p):{'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns} for p in [VAC,PHOT,SED]},
              'source_hash_scope':'Final reproducible source files; prepare and its helpers were unchanged after the successful prepare launch. Whole-file hashes are not claimed as an immutable runtime snapshot.',
              'code_sha256':{p.name:sha(p) for p in sorted(source.iterdir()) if p.suffix in ['.py','.md']},
              'validation':expected,'unit_tests':6,
              'outputs':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(OUT.iterdir()) if p.is_file()},
              'claims':'Empirical flux-completeness selection and conditional property/environment associations; no absolute completeness, causal or model-calibration claim.'}
    write_json(OUT/'MANIFEST.json',manifest)
    print('CSV tables and manifest complete',flush=True)

if __name__=='__main__':main()
