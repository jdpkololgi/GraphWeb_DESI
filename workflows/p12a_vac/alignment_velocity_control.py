"""Bounded paired velocity-scaling/RSD convention diagnostic on exposed ph000."""
import argparse
from pathlib import Path
import fitsio
import numpy as np
from mock_selection_test import digest,save
from alignment_producer_probe import ROOT
from alignment_parent_screen import BASE
from photometry_recipe_fingerprint import stats

def run(out):
    out.mkdir(parents=True,exist_ok=True)
    if (out/'VELOCITY.json').exists():raise FileExistsError(out)
    paths={b:ROOT/f'galaxy_cut_sky_{b}.fits' for b in ['N','S']}
    paths['v1']=BASE/'v1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
    a={};res={'convention_tests':{},'paired':{}}
    for name,p in paths.items():
        with fitsio.FITS(p) as f:d=f[1].read(rows=np.arange(10000),columns=['HALO_ID','CEN','RA','DEC','Z','Z_COSMO','vx','vy','vz'])
        d=d[(d['CEN']==1)&(d['Z_COSMO']>=.15)&(d['Z_COSMO']<.55)]
        keys=[(int(t['HALO_ID']),float(t['RA']),float(t['DEC'])) for t in d];assert len(keys)==len(set(keys))
        a[name]=dict(zip(keys,d))
        ra=np.deg2rad(d['RA'].astype('f8'));dec=np.deg2rad(d['DEC'].astype('f8'))
        los=d['vx']*np.cos(dec)*np.cos(ra)+d['vy']*np.cos(dec)*np.sin(ra)+d['vz']*np.sin(dec)
        apparent=(d['Z'].astype('f8')-d['Z_COSMO'])/(1+d['Z_COSMO'])*299792.458
        res['convention_tests'][name]={'rows':len(d),'cz_over_1plusz_minus_velocity_dot_sky':stats(apparent-los),
                                     'apparent_radial_velocity':stats(apparent)}
    for branch in ['N','S']:
        keys=sorted(set(a[branch])&set(a['v1']))
        old=np.array([[a['v1'][k][c] for c in ['vx','vy','vz']] for k in keys],dtype='f8')
        new=np.array([[a[branch][k][c] for c in ['vx','vy','vz']] for k in keys],dtype='f8')
        z=np.array([a['v1'][k]['Z_COSMO'] for k in keys]);assert all(a[branch][k]['Z_COSMO']==a['v1'][k]['Z_COSMO'] for k in keys)
        bands={}
        for lo,hi in [(.15,.25),(.25,.35),(.35,.45),(.45,.55)]:
            sel=(z>=lo)&(z<hi);x=old[sel];y=new[sel]
            if len(x)<2:bands[f'{lo}_{hi}']={'n':len(x)};continue
            scale=float(np.sum(x*y)/np.sum(x*x));resid=y-scale*x
            # Per-object best scalar tests the stronger necessary condition: parallel vectors.
            per=np.sum(x*y,axis=1)/np.sum(x*x,axis=1)
            bands[f'{lo}_{hi}']={'n':len(x),'least_squares_scale':scale,'residual_rms':float(np.sqrt(np.mean(resid**2))),
                'new_rms':float(np.sqrt(np.mean(y**2))),'vector_cosine_median':float(np.median(np.sum(x*y,axis=1)/(np.linalg.norm(x,axis=1)*np.linalg.norm(y,axis=1)))),
                'per_object_scalar_residual_rms':float(np.sqrt(np.mean((y-per[:,None]*x)**2)))}
        res['paired'][branch]=bands
    save(out/'VELOCITY.json',dict(**res,script_sha256=digest(__file__),sources={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in paths.values()},
       contract='First10000 raw rows per file; CEN1, .15<=Z_COSMO<.55; exact HALO_ID/RA/DEC match. Prefix-biased diagnostic, no population claim. Test v_new=s(z)*v_old, not fitting a replacement model.',
       caveat='Radial test assumes stored Cartesian velocities use the displayed sky frame and km/s; nonzero residual alone could reflect conventions. Scalar-vector residual tests do not identify underlying HOD/velocity parameters.'))
    print(res,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
