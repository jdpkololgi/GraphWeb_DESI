"""Matched occupied-sky Loa/GLAM/Holi observed-stage count screen."""
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
from flag_joint_magnitude import REAL

ROOT=Path('/global/cfs/cdirs/desi/mocks/cai/LSS/DA2/mocks')
SOURCES={name:ROOT/name/f'altmtl{num}/loa-v1/mock{num}/LSScats/BGS_BRIGHT_full_HPmapcut.dat.h5'
         for name,num in [('glam_bgs_v2',10),('holi_bgs_v2',0)]}

def run(out):
    if (out/'COMPLETE.json').exists():raise FileExistsError(out)
    common=np.load(MASK);masks={};info={}
    for name,path in SOURCES.items():
        mask=np.zeros_like(common)
        before=[path.stat().st_size,path.stat().st_mtime_ns]
        with h5py.File(path) as f:
            t=f['LSS'];n=len(t['TARGETID'])
            for start in range(0,n,500000):
                sl=slice(start,start+500000)
                mask[hp.ang2pix(256,t['RA'][sl],t['DEC'][sl],lonlat=True)]=True
            info[name]=dict(path=str(path),stat=before,rows=n,columns=list(t.keys()),
                            original_common_pixels_missing=int((common&~mask).sum()))
        masks[name]=mask
    for mask in masks.values():common &= mask
    counts={};flags={}
    for name,path in SOURCES.items():
        h=np.zeros((2,40),dtype='i8');flags[name]={}
        with h5py.File(path) as f:
            t=f['LSS'];n=len(t['TARGETID'])
            for start in range(0,n,250000):
                sl=slice(start,start+250000)
                d={k:t[k][sl] for k in ['RA','DEC','Z_not4clus','ZWARN','GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED','BGS_TARGET']}
                z=d['Z_not4clus'];c=cap(d['RA'],d['DEC'])
                keep=(z>=.15)&(z<.55)&(d['ZWARN']==0)&common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
                for i in (0,1):h[i]+=np.histogram(z[keep&(c==i)],bins=ZE)[0]
                for k in ['GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED']:
                    flags[name][k]=flags[name].get(k,0)+int((keep&~d[k].astype(bool)).sum())
                flags[name]['not_bright_bit']=flags[name].get('not_bright_bit',0)+int((keep&((d['BGS_TARGET']&2)==0)).sum())
        assert info[name]['stat']==[path.stat().st_size,path.stat().st_mtime_ns]
        counts[name]=h;print(name,h.reshape(2,4,10).sum(-1),flush=True)
    before=[REAL.stat().st_size,REAL.stat().st_mtime_ns];h=np.zeros((2,40),dtype='i8')
    with fitsio.FITS(REAL) as f:
        n=f[1].get_nrows()
        for start in range(0,n,250000):
            d=f[1].read(rows=np.arange(start,min(start+250000,n)),columns=['RA','DEC','Z_not4clus','ZWARN','DELTACHI2','SPECTYPE'])
            z=d['Z_not4clus'];c=cap(d['RA'],d['DEC'])
            keep=(z>=.15)&(z<.55)&(d['ZWARN']==0)&(d['DELTACHI2']>=25)&(np.char.strip(d['SPECTYPE'].astype('U'))=='GALAXY')&common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]
            for i in (0,1):h[i]+=np.histogram(z[keep&(c==i)],bins=ZE)[0]
        info['loa']=dict(path=str(REAL),rows=n,stat=before)
    assert before==[REAL.stat().st_size,REAL.stat().st_mtime_ns]
    counts['loa']=h
    np.save(out/'common.npy',common)
    np.savez_compressed(out/'HISTOGRAMS.npz',**counts,z_edges=ZE)
    broad={k:v.reshape(2,4,10).sum(-1) for k,v in counts.items()}
    save(out/'COMPLETE.json',dict(sources=info,counts={k:v.tolist() for k,v in broad.items()},
          mock_to_loa={k:(v/broad['loa']).tolist() for k,v in broad.items() if k!='loa'},
          flag_failures=flags,common_pixels=int(common.sum()),area_deg2=float(common.sum()*hp.nside2pixarea(256,degrees=True)),
          common_sha256=digest(out/'common.npy'),script_sha256=digest(__file__),job=os.environ.get('SLURM_JOB_ID'),
          histogram_sha256=digest(out/'HISTOGRAMS.npz'),
          notes=['One GLAM mock10 and Holi mock0 opened as development only; no independent-confirmation claim.',
                 'Observed proxy uses ZWARN0; mock DELTACHI2/SPECTYPE and photometry absent, not invented.',
                 'Intersection of occupied NSIDE256 pixels is not exact angular-mask parity.',
                 'DA3 production scripts use surveycat=DA2 and Loa; this is not a demonstrated DA3 matterhorn mock match.',
                 'No production commands executed, changes to mocks or environmental inference.']))
    print('COMPLETE internal', {k:(v/broad['loa']).tolist() for k,v in broad.items() if k!='loa'},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root)
