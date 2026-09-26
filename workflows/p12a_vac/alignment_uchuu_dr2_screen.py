"""Internal Uchuu Y3-v2.0 parent and paired clustering census; no repairs."""
import argparse
import os
from pathlib import Path
import fitsio
import h5py
import hdf5plugin
import healpy as hp
import numpy as np
from mock_selection_test import cap,digest,save
from alignment_parent_screen import MASK,ZE

ROOT=Path('/global/cfs/cdirs/desi/mocks/cai/Uchuu-SHAM/Y3-v2.0/0000')
PATHS={'paired_loa':ROOT/'BGS-BRIGHT_data/v0.1/BGS_BRIGHT_clustering.dat.fits',
       'uchuu_altmtl':ROOT/'altmtl/BGS_BRIGHT/BGS_BRIGHT_clustering.dat.h5',
       'uchuu_complete':ROOT/'complete/Uchuu-SHAM_BGS_Y3-v2.0_0000_clustering.dat.fits',
       'uchuu_any':ROOT/'complete/Uchuu-SHAM_BGS_ANY_Y3-v2.0_0000_clustering.dat.fits'}

def batches(path,columns):
    if path.suffix=='.h5':
        with h5py.File(path) as f:
            t=f['LSS'];n=len(t['RA']);cols=[k for k in columns if k in t]
            for start in range(0,n,250000):yield {k:t[k][start:start+250000] for k in cols}
    else:
        with fitsio.FITS(path) as f:
            n=f[1].get_nrows();cols=[k for k in columns if k in f[1].get_colnames()]
            for start in range(0,n,250000):
                t=f[1].read(rows=np.arange(start,min(start+250000,n)),columns=cols)
                yield {k:t[k] for k in cols}

def run(out):
    if (out/'COMPLETE.json').exists():raise FileExistsError(out)
    common=np.load(MASK);sources={}
    for key,p in PATHS.items():sources[key]=dict(path=str(p),stat=[p.stat().st_size,p.stat().st_mtime_ns])
    for name in ['paired_loa','uchuu_altmtl']:
        mask=np.zeros_like(common)
        for d in batches(PATHS[name],['RA','DEC']):mask[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]=True
        sources[name]['pixels_missing_from_initial_common']=int((np.load(MASK)&~mask).sum())
        common &= mask
    counts={};abs_hists={};weight_summaries={}
    ae=np.linspace(-27,-15,121)
    for key,p in PATHS.items():
        h=np.zeros((2,40),dtype='i8');hw=np.zeros((2,40));ah=np.zeros((2,40,120),dtype='i8');rows=0
        wmin=np.inf;wmax=-np.inf
        for d in batches(p,['RA','DEC','Z','WEIGHT','ABSMAG_R']):
            rows+=len(d['RA']);z=d['Z'];c=cap(d['RA'],d['DEC'])
            keep=(z>=.15)&(z<.55)&common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
            w=d.get('WEIGHT',np.ones(len(z)));wmin=min(wmin,float(np.min(w)));wmax=max(wmax,float(np.max(w)))
            assert np.isfinite(w).all() and np.all(w>=0)
            for i in (0,1):
                k=keep&(c==i);h[i]+=np.histogram(z[k],bins=ZE)[0]
                hw[i]+=np.histogram(z[k],bins=ZE,weights=w[k])[0]
                ah[i]+=np.histogram2d(z[k],d['ABSMAG_R'][k],bins=[ZE,ae])[0].astype('i8')
        assert sources[key]['stat']==[p.stat().st_size,p.stat().st_mtime_ns]
        sources[key]['rows']=rows;counts[key]=h;counts[key+'_weighted']=hw;abs_hists[key]=ah
        weight_summaries[key]=dict(min=wmin,max=wmax,interpretation='stored WEIGHT, never WEIGHT_FKP; missing WEIGHT uses unit count')
        print(key,h.reshape(2,4,10).sum(-1),flush=True)
    np.save(out/'common.npy',common)
    np.savez_compressed(out/'HISTOGRAMS.npz',**counts,**{k+'_zabs':v for k,v in abs_hists.items()},z_edges=ZE,abs_edges=ae)
    broad={k:v.reshape(2,4,10).sum(-1) for k,v in counts.items()}
    save(out/'COMPLETE.json',dict(sources=sources,counts={k:v.tolist() for k,v in broad.items()},
        processed_to_paired_loa=(broad['uchuu_altmtl']/broad['paired_loa']).tolist(),
        processed_to_paired_loa_weighted=(broad['uchuu_altmtl_weighted']/broad['paired_loa_weighted']).tolist(),
        weight_summaries=weight_summaries,job=os.environ.get('SLURM_JOB_ID'),script_sha256=digest(__file__),
        common_sha256=digest(out/'common.npy'),histogram_sha256=digest(out/'HISTOGRAMS.npz'),
        notes=['Newly located internal Y3-v2.0 BGS mock0000, development-only; no independent confirmation claim.',
               'Paired archived Loa v0.1 clustering selection is not the current full v2.1 DELTACHI2>=25/GALAXY sample.',
               'Weights are stored product definitions; no new n(z) reweighting. Exact assignment and spectral-success lineage remains open.',
               'Intrinsic BRIGHT/ANY counts are separate selection stages, not observed-count predictions.',
               'ABSMAG_R conventions across complete, altmtl and real data unpinned; histograms saved but not treated as physical LF parity.',
               'GALAXYID/PID exist in parents but epoch/particle/host joins and P12 estimand equivalence are unverified.']))
    print('COMPLETE DR2 Uchuu',broad['uchuu_altmtl']/broad['paired_loa'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root)
