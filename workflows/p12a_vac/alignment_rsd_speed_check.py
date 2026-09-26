"""Follow up sub-km/s radial residual with the c=300000 convention hypothesis."""
from pathlib import Path
import numpy as np
import fitsio
from alignment_producer_probe import ROOT
from alignment_parent_screen import BASE
from mock_selection_test import digest,save
from photometry_recipe_fingerprint import stats
paths={b:ROOT/f'galaxy_cut_sky_{b}.fits' for b in ['N','S']}
paths['v1']=BASE/'v1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
r={}
for name,p in paths.items():
    with fitsio.FITS(p) as f:d=f[1].read(rows=np.arange(10000),columns=['CEN','RA','DEC','Z','Z_COSMO','vx','vy','vz'])
    d=d[(d['CEN']==1)&(d['Z_COSMO']>=.15)&(d['Z_COSMO']<.55)]
    ra=np.deg2rad(d['RA'].astype('f8'));dec=np.deg2rad(d['DEC'].astype('f8'))
    v=d['vx']*np.cos(dec)*np.cos(ra)+d['vy']*np.cos(dec)*np.sin(ra)+d['vz']*np.sin(dec)
    dz=(d['Z'].astype('f8')-d['Z_COSMO'].astype('f8'))/(1+d['Z_COSMO'].astype('f8'))
    r[name]={'n':len(d),'effective_c_least_squares':float(np.sum(v*dz)/np.sum(dz*dz)),
       'residual_c299792p458':stats(dz*299792.458-v),'residual_c300000':stats(dz*300000-v)}
out=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_independent_controls_20260926/RSD_SPEED.json'
save(out,dict(results=r,script_sha256=digest(__file__),sources={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in paths.values()},contract='Same first10000 central prefix as velocity control; numerical hypothesis c=300000, no catalogue correction. Comparison tests algebraic consistency, not executed source provenance.'))
print(r)
