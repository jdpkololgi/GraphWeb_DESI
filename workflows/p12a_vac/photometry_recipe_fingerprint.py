"""Diagnostic-only numerical fingerprint of published Abacus g/r mapping.

Does not import or execute downloaded source. Independently evaluates its tables.
One uniformly spaced raw-row sample; no held-out phases or product mutations.
"""
from pathlib import Path
import hashlib, json, os
import numpy as np
import fitsio
from scipy.interpolate import splrep, splev
from astropy.cosmology import FlatLambdaCDM
OUT=Path(__file__).resolve().parents[2]/'docs/evidence/p12a_photometry_recipe_20260926'
RAW=Path('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v0.1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph006.fits')
def k(z,c,band,cubic=True,extend=True):
    t=np.loadtxt(OUT/'source'/f'k_corr_{band}band_z01.dat'); cm=t[:,7]; cc=np.clip(c,cm.min(),cm.max())
    def poly(zz):
        co=[splev(cc,splrep(cm,t[:,i])) if cubic else np.interp(cc,cm,t[:,i]) for i in range(2,6)]
        return sum(a*(zz-.1)**p for a,p in zip(co,[4,3,2,1]))+t[0,6]
    val=poly(z)
    if extend:
        # Source uses linear interpolation of secant coefficients in colour,
        # even when the polynomial branch uses cubic interpolation.
        def knotpoly(zz): return sum(t[:,i]*(zz-.1)**p for i,p in zip(range(2,6),[4,3,2,1]))+t[0,6]
        slope=(knotpoly(.5)-knotpoly(.4))/.1
        line=np.interp(cc,cm,slope)*z+np.interp(cc,cm,knotpoly(.4)-slope*.4)
        val=np.where(z>.5,line,val)
    return val

def stats(a):
    a=a[np.isfinite(a)]
    return dict(n=len(a),median=float(np.median(a)),p01_p99=np.quantile(a,[.01,.99]).tolist(),max_abs=float(np.max(np.abs(a))),rms=float(np.sqrt(np.mean(a*a)))) if len(a) else dict(n=0)
def main():
    with fitsio.FITS(RAW) as f:
        n=f[1].get_nrows(); rows=np.arange(0,n,320)
        d=f[1].read(rows=rows,columns=['Z','Z_COSMO','R_MAG_APP','R_MAG_ABS','G_R_REST','G_R_OBS'])
    np.savez_compressed(OUT/'RAW_SAMPLE.npz',data=d,rows=rows)
    z=d['Z']; c=d['G_R_REST']; keep=(z>=.15)&(z<.55)&(d['R_MAG_APP']<19.5)
    out=dict(raw=str(RAW),raw_rows=n,sample_rows=len(d),stride=320,job=os.environ.get('SLURM_JOB_ID'),source_url='https://github.com/amjsmith/shared_code',source_commit='a40a5d3be570a6ce9d5a0698268b04864d4d3af2',variants={})
    # Native-like distance approximation; colour test is distance independent.
    cosmo=FlatLambdaCDM(H0=67.36,Om0=(.02237+.12+.0006442)/.6736**2,Tcmb0=0)
    grid=np.linspace(0,1,20001); dl=cosmo.luminosity_distance(grid).value*cosmo.h
    for zn in ['Z','Z_COSMO']:
        zz=d[zn]
        for cubic in [True,False]:
            kr=k(zz,c,'r',cubic);kg=k(zz,c,'g',cubic)
            dc=d['G_R_OBS']-(c+kg-kr)
            dm=d['R_MAG_APP']-(d['R_MAG_ABS']+5*np.log10(np.maximum(np.interp(zz,grid,dl),1e-20))+25+kr)
            out['variants'][f'{zn}_cubic{cubic}']={f'{lo:.2f}_{hi:.2f}':dict(colour_residual=stats(dc[keep&(z>=lo)&(z<hi)]),r_residual_approx_distance=stats(dm[keep&(z>=lo)&(z<hi)])) for lo,hi in [(.15,.25),(.25,.35),(.35,.45),(.45,.5),(.5,.55)]}
    kr=k(z,c,'r'); residual=d['R_MAG_APP']-(d['R_MAG_ABS']+5*np.log10(np.maximum(np.interp(z,grid,dl),1e-20))+25+kr)
    fit=np.polyfit(z[keep]-.1,residual[keep],1)
    out['r_residual_linear_fit']=dict(slope=float(fit[0]),intercept=float(fit[1]),after_fit=stats((residual-np.polyval(fit,z-.1))[keep]),after_Q07=stats((residual+.7*(z-.1))[keep]),after_Q08=stats((residual+.8*(z-.1))[keep]))
    out['distance_caveat']='Flat LCDM Om=(omega_b+omega_cdm+omega_ncdm)/h^2, radiation neglected; r residual not an exact cosmology replay. Colour identity is distance independent.'
    out['highz_branch_diagnostic']=dict(r_linear_minus_polynomial=stats((k(z,c,'r')-k(z,c,'r',extend=False))[keep&(z>.5)]),g_linear_minus_polynomial=stats((k(z,c,'g')-k(z,c,'g',extend=False))[keep&(z>.5)]))
    out['hashes']={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((OUT/'source').iterdir())}
    out['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'FINGERPRINT.json').write_text(json.dumps(out,indent=2)+'\n'); print(json.dumps(out,indent=2))
if __name__=='__main__':main()
