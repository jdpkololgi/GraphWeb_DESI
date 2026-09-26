"""Chunked, restartable intrinsic-parent version screen; no catalogue repairs."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
import fitsio
import healpy as hp
from mock_selection_test import cap, digest, save

BASE = Path('/global/cfs/cdirs/desi/cosmosim/SecondGenMocks/AbacusSummit/CutSky/BGS')
MASK = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/p12a_mock_test_ph006_20260926_v1/common.npy')
ZE = np.linspace(.15, .55, 41)
ME = np.linspace(12, 20.2, 83)
CE = np.linspace(-1, 4, 101)
AE = np.linspace(-26, -14, 121)


def run(root, version, chunk):
    source = BASE / version / 'z0.200/cutsky_BGS_z0.200_AbacusSummit_base_c000_ph000.fits'
    receipt = root / f'{version}.json'
    if receipt.exists():
        old = json.loads(receipt.read_text())
        assert old['source_stat'] == [source.stat().st_size, source.stat().st_mtime_ns]
        assert old['script_sha256'] == digest(__file__)
        assert old['histogram_sha256'] == digest(root / f'{version}.npz')
        print('already complete', version, flush=True)
        return
    common = np.load(MASK)
    before = [source.stat().st_size, source.stat().st_mtime_ns]
    hist = dict(z=np.zeros((2,40), dtype='i8'),
                zr=np.zeros((2,40,82), dtype='i8'),
                zcolour=np.zeros((2,40,100), dtype='i8'),
                zabs=np.zeros((2,40,120), dtype='i8'),
                central=np.zeros((2,40), dtype='i8'))
    parent = np.zeros((2,40), dtype='i8')
    invalid = 0
    with fitsio.FITS(source) as f:
        n = f[1].get_nrows()
        columns = f[1].get_colnames()
        for start in range(0,n,chunk):
            d = f[1].read(rows=np.arange(start,min(start+chunk,n)),
                          columns=['RA','DEC','Z','R_MAG_APP','R_MAG_ABS','G_R_OBS','CEN'])
            finite = np.isfinite(d['RA']) & np.isfinite(d['DEC']) & np.isfinite(d['Z']) & np.isfinite(d['R_MAG_APP'])
            invalid += int((~finite).sum())
            d = d[finite & (d['Z']>=ZE[0]) & (d['Z']<ZE[-1])]
            d = d[common[hp.ang2pix(256,d['RA'],d['DEC'],lonlat=True)]]
            c = cap(d['RA'],d['DEC'])
            for i in (0,1):
                parent[i] += np.histogram(d['Z'][c==i], bins=ZE)[0]
                t = d[(c==i) & (d['R_MAG_APP']>=12) & (d['R_MAG_APP']<19.5)]
                hist['z'][i] += np.histogram(t['Z'], bins=ZE)[0]
                hist['central'][i] += np.histogram(t['Z'][t['CEN']==1], bins=ZE)[0]
                for name,col,edges in [('zr','R_MAG_APP',ME),('zcolour','G_R_OBS',CE),('zabs','R_MAG_ABS',AE)]:
                    hist[name][i] += np.histogram2d(t['Z'],t[col],bins=[ZE,edges])[0].astype('i8')
            if start % (chunk*8)==0:
                print(version,start,n,flush=True)
    assert before == [source.stat().st_size,source.stat().st_mtime_ns]
    assert np.array_equal(hist['zr'].sum(-1),hist['z'])
    assert np.all(hist['central'] <= hist['z'])
    np.savez_compressed(root/f'{version}.npz',**hist,parent=parent,z_edges=ZE,r_edges=ME,colour_edges=CE,abs_edges=AE)
    save(receipt,dict(source=str(source),source_stat=before,total_rows=n,columns=columns,
        selected=int(hist['z'].sum()),invalid_coordinates_z_r=invalid,
        counts=hist['z'].reshape(2,4,10).sum(-1).tolist(),
        central_counts=hist['central'].reshape(2,4,10).sum(-1).tolist(),
        job=os.environ.get('SLURM_JOB_ID'),script_sha256=digest(__file__),
        histogram_sha256=digest(root/f'{version}.npz'),common_mask_sha256=digest(MASK),
        contract={'phase':'ph000','caps':['SGC','NGC'],'z':'0.15<=Z<0.55 RSD',
                  'magnitude':'12<=R_MAG_APP<19.5','IN_Y_flags':'none in either version',
                  'population':'intrinsic raw parent, not observed Loa',
                  'photometry_equivalent_to_Legacy':False,'replicate_independence':False},
        validation={'source_unchanged':True,'magnitude_histogram_conserves_counts':True}))
    print('COMPLETE',version,hist['z'].reshape(2,4,10).sum(-1),flush=True)


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--root',type=Path,required=True)
    p.add_argument('--version',choices=['v0.1','v1'],required=True)
    p.add_argument('--chunk',type=int,default=500000)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True)
    run(a.root,a.version,a.chunk)
