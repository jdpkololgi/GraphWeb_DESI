"""Recompute published primary curves without calling the estimator."""
import json
import os
import socket
import sys
from pathlib import Path
import numpy as np

ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_mass_complete_20261016_v1')


def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Requires compute node')
    full='--full' in sys.argv
    if full:
        import fitsio
        a=fitsio.read(ROOT/'LOA_P12A_MASS_COMPLETE_zdependent.fits')
    else:a=np.load(ROOT/'analysis_input.npy',allow_pickle=False)
    result=json.loads((ROOT/('relationships_full.json' if full else 'relationships.json')).read_text())
    checks=0
    for width in [.1,.05]:
        for threshold,key in [(-11.,'quenched_ssfr11'),(-10.5,'quenched_ssfr10p5'),(-11.5,'quenched_ssfr11p5')]:
            valid=(np.ones(len(a),bool) if full else a['SELECT_FIXED020030']) & a['VALID_SFR']
            d=a[valid];m=np.log10(d['MASS_CG_15'])
            # Avoid logarithmic zero handling entirely in this reconstruction.
            y=(d['SFR_CG_15']/d['MASS_CG_15'] < 10**threshold).astype(float)
            cells=np.floor((m-6)/width).astype(int)*100+np.floor((d['Z']-.2+1e-9)/.01).astype(int)
            rec=result['zdependent_full' if full else 'fixed'][str(width)][key]
            assert rec['status']=='measured'
            ids,inv=np.unique(cells,return_inverse=True)
            keep=np.isin(ids,rec['common_cell_ids']);ref=np.bincount(inv)[keep].astype(float);ref/=ref.sum()
            np.testing.assert_allclose(ref,rec['reference_weights'],rtol=0,atol=1e-14)
            assert int(keep[inv].sum())==rec['rows']
            p=d['P_WEB'].astype(float)
            for name,w in [('soft_web',p),('hard_web',np.eye(4)[np.argmax(p,axis=1)])]:
                mu=[];raw=[]
                for j in range(4):
                    count=np.bincount(inv,weights=w[:,j],minlength=len(ids))[keep]
                    numerator=np.bincount(inv,weights=w[:,j]*y,minlength=len(ids))[keep]
                    mu.append(np.dot(ref,numerator/count));raw.append(numerator.sum()/count.sum())
                published=rec['methods'][name]
                np.testing.assert_allclose(mu,published['mean'],rtol=0,atol=2e-10)
                np.testing.assert_allclose(raw,published['raw_on_common_rows'],rtol=0,atol=2e-10)
                assert abs((mu[-1]-mu[0])-published['contrast'])<2e-10
                ci=np.array(published['ci95']);assert ci.shape==(2,4) and np.isfinite(ci).all() and (ci[0]<=ci[1]).all()
                checks+=1
    with (ROOT/('relationship_validation_full.json' if full else 'relationship_validation.json')).open('x') as f:
        json.dump({'passed':checks,'scope':'Independent binned weighted means, same-row raw means and contrasts for three quenching thresholds and two mass bin widths, both class estimators.'},f,indent=2)
    print('Independent curve reconstructions passed:',checks,flush=True)

if __name__=='__main__':main()
