"""Independent reconstruction using stored support and direct weighted histograms."""
import os
import socket
import json
import numpy as np
import fitsio
from run import OUT, WINDOWS, EDGES, sha, write

def main():
    assert os.environ.get('SLURM_JOB_ID') and socket.gethostname().startswith('nid')
    manifest=json.loads((OUT/'MANIFEST.json').read_text())
    for name,h in manifest['files'].items():assert sha(OUT/name)==h
    results=[]
    for key,*_ in WINDOWS:
        r=json.loads((OUT/f'{key}.json').read_text())
        a=fitsio.read(r['source'],columns=['MASS_CG_15','Z','GR_OBS','P_WEB'])
        m=np.log10(a['MASS_CG_15']);z=a['Z'];y=a['GR_OBS'];p=a['P_WEB'].astype(float)
        for grid,g in r['grids'].items():
            cell=np.column_stack((np.floor(m/g['dm']+1e-9).astype(int),np.floor(z/g['dz']+1e-9).astype(int)))
            lookup={tuple(v):i for i,v in enumerate(g['cell_ids'])}
            c=np.fromiter((lookup.get(tuple(v),-1) for v in cell),dtype=int,count=len(cell));keep=c>=0
            target=np.array(g['reference_weights']);nc=len(target)
            assert keep.sum()==g['rows']
            assert np.allclose(np.bincount(c[keep],minlength=nc)/keep.sum(),target)
            for method,w in [('hard',np.eye(4)[p.argmax(1)]),('soft',p)]:
                rec=g['methods'][method];actual=[];rawmeans=[];rawcommon=[];mm=[];zz=[]
                for j in range(4):
                    rawmeans.append(float(np.average(y,weights=w[:,j])))
                    rawcommon.append(float(np.average(y[keep],weights=w[keep,j])))
                    den=np.bincount(c[keep],weights=w[keep,j],minlength=nc)
                    weight=w[keep,j]*target[c[keep]]/den[c[keep]]
                    hist=np.histogram(y[keep],bins=np.r_[-np.inf,EDGES,np.inf],weights=weight)[0]
                    assert np.allclose(hist,rec['matched'][j],atol=1e-10)
                    actual.append(float(np.average(y[keep],weights=weight)))
                    mm.append(float(np.average(m[keep],weights=weight)));zz.append(float(np.average(z[keep],weights=weight)))
                    assert np.isclose(mm[-1],rec['balance'][j]['mean_logmass'])
                    assert np.isclose(zz[-1],rec['balance'][j]['mean_z'])
                assert np.allclose(actual,rec['mean_colour'],atol=1e-12)
                assert np.isclose(actual[3]-actual[0],rec['cluster_minus_void'])
                assert rec['bootstrap_valid']==256
                results.append(dict(window=key,grid=grid,method=method,raw_all_mean_colour=rawmeans,
                    raw_common_mean_colour=rawcommon,matched_mean_colour=actual,
                    matched_mass_mean_range=float(np.ptp(mm)),matched_z_mean_range=float(np.ptp(zz)),
                    max_cell_balance_error=max(b['max_cell_error'] for b in rec['balance']),
                    neff=[b['neff'] for b in rec['balance']]))
    write(OUT/'INDEPENDENT_VALIDATION.json',dict(histograms_checked=48,all_bootstraps_valid=True,results=results))
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
