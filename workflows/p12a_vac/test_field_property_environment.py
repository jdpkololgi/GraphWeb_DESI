"""Synthetic end-to-end figure smoke, never observational evidence."""
import argparse
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
import numpy as np
from workflows.p12a_vac import field_property_environment as analysis


def main():
    if 'SLURM_JOB_ID' not in os.environ:raise RuntimeError('Run figure smoke on compute')
    analysis.tests()
    rng=np.random.default_rng(781);n=24000
    with tempfile.TemporaryDirectory(prefix='loa-property-fixture-') as temp:
        root=Path(temp);out=root/'out';z=rng.uniform(.15,.55,n);m=rng.uniform(9.5,11.5,n)
        labs=rng.integers(0,4,size=(8,n),dtype='u1');p=np.stack([(labs==j).mean(0) for j in range(4)],1)
        mass=10**m;sfr=mass*10**rng.normal(-10.7,.7,n)
        np.savez_compressed(root/'vac_properties.npz',TARGETID=np.arange(n),Z=z,RA=rng.uniform(0,360,n),DEC=rng.uniform(-20,60,n),
            CAP=rng.integers(0,2,n),P_WEB=np.full((n,4),.25),SUPPORTED=np.ones(n,bool),CIGALE_MATCHED=np.ones(n,bool),
            CIGALE_SPECTYPE=np.full(n,'GALAXY'),MASS_CG_15=mass,SFR_CG_15=sfr,MASS_CG_5=mass*1.05,SFR_CG_5=sfr*.98)
        np.savez_compressed(root/'FIELD_MARGINALS.npz',p_web=p,p_web_threshold0=p[:,::-1],classes02=labs,field_eligible=np.ones(n,bool))
        (root/'PLAN.json').write_text(json.dumps({'properties_sha256':analysis.digest(root/'vac_properties.npz'),'property_product':{'fixture':True}}))
        (root/'ATLAS_COMPLETE.json').write_text(json.dumps({'passed':True,'rows':n,'draws':8,'output':{'path':str(root/'FIELD_MARGINALS.npz'),'sha256':analysis.digest(root/'FIELD_MARGINALS.npz')}}))
        analysis.main(SimpleNamespace(root=root,out=out))
        r=json.loads((out/'RESULTS.json').read_text())
        assert len(r['figures'])==8 and len(list(out.glob('*.png')))==8
        assert r['rows']==n and r['common_both_fits']==n
        assert set(r['threshold_sensitivity'])=={'logssfr','low_ssfr'}
        print('Synthetic full figure pipeline passed; fixture outputs removed.',flush=True)


if __name__=='__main__':main()
