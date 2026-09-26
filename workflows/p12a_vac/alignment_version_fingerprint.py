"""Compare numerical K/E fingerprints across the exposed ph000 versions."""
import argparse
import os
from pathlib import Path
import fitsio
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from photometry_recipe_fingerprint import k,stats,OUT
from alignment_parent_screen import BASE
from mock_selection_test import digest,save

def run(root):
    cosmo=FlatLambdaCDM(H0=67.36,Om0=(.02237+.12+.0006442)/.6736**2,Tcmb0=0)
    grid=np.linspace(0,1,20001);dl=cosmo.luminosity_distance(grid).value*cosmo.h
    result={}
    for version in ['v0.1','v1']:
        p=BASE/version/'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
        with fitsio.FITS(p) as f:
            rows=np.arange(0,f[1].get_nrows(),320)
            d=f[1].read(rows=rows,columns=['Z','Z_COSMO','G_R_REST','G_R_OBS','R_MAG_APP','R_MAG_ABS'])
        z=d['Z'];c=d['G_R_REST'];keep=(z>=.15)&(z<.55)&(d['R_MAG_APP']<19.5)
        values={}
        for name in ['Z','Z_COSMO']:
            zz=d[name];kr=k(zz,c,'r');kg=k(zz,c,'g')
            residual=d['R_MAG_APP']-d['R_MAG_ABS']-5*np.log10(np.maximum(np.interp(zz,grid,dl),1e-20))-25-kr
            fit=np.polyfit(zz[keep]-.1,residual[keep],1)
            values[name]=dict(fit_slope_intercept=fit.tolist(),r_after_fit=stats((residual-np.polyval(fit,zz-.1))[keep]),
                             colour_residual=stats((d['G_R_OBS']-c-kg+kr)[keep]))
        result[version]=dict(source=str(p),sample_rows=len(d),selected_rows=int(keep.sum()),values=values)
    save(root/'FINGERPRINT.json',dict(results=result,job=os.environ.get('SLURM_JOB_ID'),
         script_sha256=digest(__file__),evaluator_sha256=digest(Path(__file__).with_name('photometry_recipe_fingerprint.py')),
         table_sha256={p.name:digest(p) for p in (OUT/'source').glob('k_corr*')},
         caveat='Native-like distance approximation; numerical identity is not proof of executed LF parameters or production commit.'))
    print(result,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root)
