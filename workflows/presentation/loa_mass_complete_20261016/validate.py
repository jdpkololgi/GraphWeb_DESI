"""Independent source-row and selection checks for the generated sub-VACs."""
import os
import socket
from pathlib import Path
import json
import numpy as np
import fitsio

ROOT=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_mass_complete_20261016_v1')
SOURCE=Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_loa_cigale_20260929_v1/DESI_LOA_P12A_CIGALE_FULL_VAC.fits')


def main():
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('Validation requires a compute node')
    curve=json.loads((ROOT/'mass_limits.json').read_text())
    summary=json.loads((ROOT/'selection_summary.json').read_text())
    parent=np.load(ROOT/'bright_parent.npz')
    checks=[]
    # Independently reconstruct conservative limits, sorting the rescaled
    # faint population and explicitly interpolating the percentile order stats.
    for b in curve:
        estimates=[]
        base=(parent['z']>=b['zlo'])&(parent['z']<b['zhi'])
        sf=parent['sfr'];m=parent['mass'];off=parent['r']-parent['rlim']
        q=(sf>=0)&np.isfinite(sf)&(sf/10**m<1e-11)
        for name,sel in [('all',base),('quiescent',base&q),('starforming',base&~q&(sf>=0)&np.isfinite(sf))]:
            ref=b['populations'][name]
            assert int(sel.sum())==ref['rows']
            if ref['usable']:
                threshold=np.percentile(off[sel],80)
                v=np.sort((m+.4*off)[sel&(off>=threshold)])
                t=(len(v)-1)*.95;lo=int(np.floor(t));hi=int(np.ceil(t))
                est=v[lo]+(v[hi]-v[lo])*(t-lo)
                assert abs(est-ref['p95'])<1e-10
                estimates.append(est)
                checks.append('limiting_mass_'+name+'_'+str(b['zlo']))
        if b['usable']:assert abs(max(estimates)-b['logmass_limit'])<1e-10
    for name in ['zdependent','fixed020030']:
        path=ROOT/f'LOA_P12A_MASS_COMPLETE_{name}.fits'
        a=fitsio.read(path);assert len(np.unique(a['TARGETID']))==len(a)
        assert a['SUPPORTED'].all() and a['CIGALE_MATCHED'].all() and not a['CIGALE_AMBIGUOUS'].any()
        assert np.all(a['R_DERED']<a['R_LIMIT'])
        mass=np.log10(a['MASS_CG_15'])
        if name=='fixed020030':
            assert ((a['Z']>=.2)&(a['Z']<.3)&(mass>=summary['fixed_logmass'])).all()
            assert len(a)==summary['fixed020030_rows']
        else:
            assert (mass>=a['LOGM_LIMIT95']).all()
            assert len(a)==summary['zdependent_rows']
        k=np.linspace(0,len(a)-1,min(257,len(a)),dtype=int)
        rows=a['SOURCE_ROW'][k];order=np.argsort(rows)
        source=fitsio.read(SOURCE,rows=rows[order])[np.argsort(order)]
        for col in source.dtype.names:
            actual=a[col][k]
            if source[col].dtype.kind in 'fc':np.testing.assert_allclose(actual,source[col],rtol=0,atol=0,equal_nan=True)
            else:np.testing.assert_array_equal(actual,source[col])
        checks.append(name+'_all_original_columns_257_source_rows')
        sfr=a['SFR_CG_15'];valid=np.isfinite(sfr)&(sfr>=0)
        assert np.array_equal(valid,a['VALID_SFR'])
        assert np.array_equal(a['QUENCHED_SSFR11'][valid],(sfr[valid]/a['MASS_CG_15'][valid]<1e-11).astype(int))
        assert (a['QUENCHED_SSFR11'][~valid]==-1).all()
        checks.append(name+'_selection_unique_ids_flux_quenching')
    result={'passed':len(checks),'checks':checks,'source_row_checks_per_catalogue':257,
            'scope':'Selection arithmetic, all original columns at source rows, unique IDs and quenching flags; not a validation of absolute survey completeness or environment predictions.'}
    with (ROOT/'validation.json').open('x') as f:json.dump(result,f,indent=2)
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
