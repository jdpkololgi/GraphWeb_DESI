"""Angular erosion and tile-count controls on original observed mock versus Loa."""
import argparse,os,json
from pathlib import Path
import fitsio,h5py,hdf5plugin
import healpy as hp
import numpy as np
from mock_selection_test import cap,digest,save
from alignment_parent_screen import MASK,ZE
from flag_joint_magnitude import REAL,BASE
MOCK=BASE/'altmtl6/loa-v1/mock6/LSScats/BGS_BRIGHT_full_HPmapcut.dat.h5'

def erode(m):
    result=np.zeros_like(m);ix=np.flatnonzero(m)
    nei=hp.get_all_neighbours(256,ix,nest=True)
    result[ix]=np.all((nei<0)|m[np.maximum(nei,0)],axis=0)
    return result

def accumulate(hist,d,masks,real):
    z=d['Z_not4clus'];ok=(d['ZWARN']==0)&(z>=.15)&(z<.55)
    if real:ok&=(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')
    ix=np.flatnonzero(ok);z=z[ix];ra=d['RA'][ix];dec=d['DEC'][ix]
    pix=hp.ang2pix(512,ra,dec,lonlat=True,nest=True);c=cap(ra,dec)
    tile=np.clip(d['NTILE'][ix],1,3).astype(int)-1
    for name,m in masks.items():
        k=m[pix]
        for j in (0,1):
            for nt in range(3):hist[name][j,nt]+=np.histogram(z[k&(c==j)&(tile==nt)],ZE)[0]

def run(root):
    root.mkdir(parents=True,exist_ok=True)
    if (root/'COMPLETE.json').exists():raise FileExistsError(root)
    before={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [REAL,MOCK]}
    orig=hp.reorder(np.load(MASK),r2n=True).astype(bool);one=erode(orig);two=erode(one)
    masks={'original':np.repeat(orig,4),'erode1':np.repeat(one,4),'erode2':np.repeat(two,4)}
    hist={};sources={}
    cols=['RA','DEC','Z_not4clus','ZWARN','NTILE']
    for key,p in [('loa',REAL),('mock',MOCK)]:
        h={name:np.zeros((2,3,40),dtype='i8') for name in masks}
        if key=='loa':
            with fitsio.FITS(p) as f:
                n=f[1].get_nrows()
                for start in range(0,n,250000):accumulate(h,f[1].read(rows=np.arange(start,min(start+250000,n)),columns=cols+['DELTACHI2','SPECTYPE']),masks,True)
        else:
            with h5py.File(p) as f:
                t=f['LSS'];n=len(t['RA'])
                for start in range(0,n,250000):accumulate(h,{k:t[k][start:start+250000] for k in cols},masks,False)
        hist.update({key+'_'+k:v for k,v in h.items()});print(key,'done',flush=True)
    baseline=np.load(Path(__file__).resolve().parents[2]/'docs/evidence/p12a_joint_flags_20260926/HISTOGRAMS.npz')
    assert np.array_equal(hist['loa_original'].sum(1),baseline['combinations'][...,0])
    assert np.array_equal(hist['mock_original'].sum(1),baseline['mock_observed'].sum(-1))
    for key in ['loa','mock']:
        assert np.all(hist[key+'_erode2']<=hist[key+'_erode1'])
        assert np.all(hist[key+'_erode1']<=hist[key+'_original'])
    broad={k:v.sum(1).reshape(2,4,10).sum(-1) for k,v in hist.items()}
    assert before=={str(p):[p.stat().st_size,p.stat().st_mtime_ns] for p in [REAL,MOCK]}
    np.savez_compressed(root/'HISTOGRAMS.npz',**hist,z_edges=ZE)
    save(root/'COMPLETE.json',dict(counts={k:v.tolist() for k,v in broad.items()},loa_to_mock={k:(broad['loa_'+k]/broad['mock_'+k]).tolist() for k in masks},
        mask_pixel_counts={k:int(v.sum()) for k,v in masks.items()},sources=before,job=os.environ.get('SLURM_JOB_ID'),
        script_sha256=digest(__file__),helper_sha256={n:digest(Path(__file__).with_name(n)) for n in ['mock_selection_test.py','alignment_parent_screen.py','flag_joint_magnitude.py']},
        histogram_sha256=digest(root/'HISTOGRAMS.npz'),mask_sha256=digest(MASK),
        validation={'baseline_reproduced':True,'erosion_nested':True},
        contract='Original exposed ph006 observed-stage altMTL Loa mock, ZWARN0, versus full Loa quality25+GALAXY, .15<=z<.55. No weights, no new phases. NTILE <=1,2,>=3 strata. Erode one/two rings of original NSIDE256 occupied support.',
        caveat='Erosion is a sensitivity control, not exact matched masks or an observation model. Tile-count strata do not match density, competition or subpixel selection.'))
    print('COMPLETE boundary', {k:(broad['loa_'+k]/broad['mock_'+k]).tolist() for k in masks},flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);run(p.parse_args().root)
