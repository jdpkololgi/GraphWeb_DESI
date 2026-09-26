"""Independent producer N/S ph000 censuses; never concatenate unverified branches."""
import argparse
import os
from pathlib import Path
import fitsio
import healpy as hp
import numpy as np
from alignment_producer_probe import ROOT
from alignment_parent_screen import MASK, ZE, ME, CE
from mock_selection_test import cap, digest, save

def run(out):
    out.mkdir(parents=True, exist_ok=True)
    common = np.load(MASK)
    for p in sorted(ROOT.glob('*.fits')) + sorted((ROOT/'forFA').glob('*.fits')):
        receipt=out/(p.stem+'.json'); product=out/(p.stem+'.npz')
        if receipt.exists():
            raise FileExistsError(receipt)
        before=[p.stat().st_size,p.stat().st_mtime_ns]
        counts=np.zeros((2,40),dtype='i8'); zr=np.zeros((2,40,82),dtype='i8'); zc=np.zeros((2,40,100),dtype='i8')
        flags={}; invalid=0
        with fitsio.FITS(p) as f:
            n=f[1].get_nrows(); names=f[1].get_colnames(); zn='Z' if 'Z' in names else 'RSDZ'
            cols=['RA','DEC',zn,'R_MAG_APP','G_R_OBS']+(['BGS_TARGET'] if 'BGS_TARGET' in names else [])
            for start in range(0,n,500000):
                d=f[1].read(rows=np.arange(start,min(start+500000,n)),columns=cols)
                if 'BGS_TARGET' in names:
                    u,c=np.unique(d['BGS_TARGET'],return_counts=True)
                    for a,b in zip(u,c):flags[str(int(a))]=flags.get(str(int(a)),0)+int(b)
                ok=np.isfinite(d['RA'])&np.isfinite(d['DEC'])&np.isfinite(d[zn])&np.isfinite(d['R_MAG_APP'])
                invalid+=int((~ok).sum())
                d=d[ok&(d[zn]>=.15)&(d[zn]<.55)&(d['R_MAG_APP']>=12)&(d['R_MAG_APP']<19.5)]
                d=d[common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]]
                hemis=cap(d['RA'],d['DEC'])
                for i in [0,1]:
                    t=d[hemis==i]
                    counts[i]+=np.histogram(t[zn],ZE)[0]
                    zr[i]+=np.histogram2d(t[zn],t['R_MAG_APP'],[ZE,ME])[0].astype('i8')
                    zc[i]+=np.histogram2d(t[zn],t['G_R_OBS'],[ZE,CE])[0].astype('i8')
                if start%5000000==0:print(p.stem,start,n,flush=True)
        assert before==[p.stat().st_size,p.stat().st_mtime_ns]
        assert np.array_equal(zr.sum(-1),counts)
        np.savez_compressed(product,z=counts,zr=zr,zcolour=zc,z_edges=ZE,r_edges=ME,colour_edges=CE)
        save(receipt,dict(path=str(p),source_stat=before,total_rows=n,invalid=invalid,counts=counts.reshape(2,4,10).sum(-1).tolist(),
             flags_all_rows=flags,script_sha256=digest(__file__),histogram_sha256=digest(product),mask_sha256=digest(MASK),
             helper_sha256={s:digest(Path(__file__).with_name(s)) for s in ['mock_selection_test.py','alignment_parent_screen.py','alignment_producer_probe.py']},
             job=os.environ.get('SLURM_JOB_ID'),contract='Separate N/S products, not additive; common occupied NSIDE256 sky, galactic SGC/NGC, 12<=r<19.5, .15<=RSD z<.55. Numerical magnitude cut only; no claim of Legacy parity or observed quality.'))
        print('COMPLETE',p.stem,counts.reshape(2,4,10).sum(-1).tolist(),flags,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    run(p.parse_args().root)
