"""Bounded deterministic ph000 photometry fingerprint; no catalogue mutation."""
import argparse, json, sys, hashlib
from pathlib import Path
import numpy as np
import fitsio
from scipy.integrate import cumulative_trapezoid

def stats(x):
    x=np.asarray(x,dtype=float)
    return dict(n=int(x.size), finite=int(np.isfinite(x).sum()), rms=float(np.sqrt(np.mean(x*x))), median=float(np.median(x)), maxabs=float(np.max(np.abs(x))))

def main():
    p=argparse.ArgumentParser(); p.add_argument('--source',required=True); p.add_argument('--output',required=True); a=p.parse_args()
    sys.path.insert(0,a.source)
    from hodpy.k_correction import DESI_KCorrection, DESI_KCorrection_color
    path='/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS/v1/z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
    with fitsio.FITS(path) as f:
        n=f[1].get_nrows(); rows=np.arange(0,n,1000); d=f[1].read(rows=rows)
    d=d[(d['Z']>=.15)&(d['Z']<.55)]
    z=d['Z'].astype(float); zc=d['Z_COSMO'].astype(float); c=d['G_R_REST'].astype(float)
    # Approximate native c000 flat LCDM; isolate exact colour from distance-dependent r.
    grid=np.linspace(0,1,100001); om=(.02237+.1200+.00064420)/.6736**2
    chi=cumulative_trapezoid(299792.458/100/np.sqrt(om*(1+grid)**3+1-om),grid,initial=0)
    result={'path':path,'total_rows':n,'stride':1000,'selected_rows':len(d),'distance':'approximate flat LCDM c000, Mpc/h; not exact CLASS replay','omega_m':om,'source_hashes':{str(x.relative_to(a.source)):hashlib.sha256(x.read_bytes()).hexdigest() for x in Path(a.source).rglob('*') if x.is_file() and '__pycache__' not in str(x)},'variants':{}}
    residuals={}
    for ps in ['N','S']:
        k=DESI_KCorrection('r',ps); kc=DESI_KCorrection_color(ps)
        col=kc.observer_frame_colour(z,np.clip(c,-3.9,3.9))-d['G_R_OBS']
        for convention,dl in [('zobs',(1+z)*np.interp(z,grid,chi)),('mixed',(1+z)**2/(1+zc)*np.interp(zc,grid,chi))]:
            r=d['R_MAG_ABS']+5*np.log10(dl)+25+k.k(z,c)-.67*(z-.1)-d['R_MAG_APP']
            residuals[ps+'_'+convention]=(col,r)
            result['variants'][ps+'_'+convention]={'colour':stats(col),'r':stats(r),'shells':[{'zlo':lo,'zhi':hi,'colour':stats(col[(z>=lo)&(z<hi)]),'r':stats(r[(z>=lo)&(z<hi)])} for lo,hi in zip([.15,.25,.35,.45],[.25,.35,.45,.55])]}
    cn,rn=residuals['N_zobs']; cs,rs=residuals['S_zobs']
    north=np.abs(cn)<np.abs(cs)
    result['inferred_photsys']={'note':'chosen by minimum colour residual; diagnostic inference, not known provenance', 'north':int(north.sum()), 'south':int((~north).sum()), 'colour':stats(np.where(north,cn,cs)), 'r':stats(np.where(north,rn,rs)), 'north_dec':np.percentile(d['DEC'][north],[0,50,100]).tolist(), 'south_dec':np.percentile(d['DEC'][~north],[0,50,100]).tolist(), 'bright_colour':stats(np.where(north,cn,cs)[d['R_MAG_APP']<19.5]), 'bright_r':stats(np.where(north,rn,rs)[d['R_MAG_APP']<19.5])}
    fixed=d['DEC']>32.375
    result['fixed_dec_split']={'threshold':32.375,'disagreements_with_colour_inference':int(np.sum(fixed!=north)), 'colour':stats(np.where(fixed,cn,cs)), 'r_zobs':stats(np.where(fixed,rn,rs)), 'r_mixed':stats(np.where(fixed,residuals['N_mixed'][1],residuals['S_mixed'][1]))}
    assert result['fixed_dec_split']['colour']['finite']==len(d)
    assert sum(v['colour']['n'] for v in result['variants']['S_zobs']['shells'])==len(d)
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:{s:v[s] for s in ['colour','r']} for k,v in result['variants'].items()},indent=2))
if __name__=='__main__': main()
