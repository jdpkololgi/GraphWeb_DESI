"""Uchuu SV3 ensemble census against the paper's archived SV3 n(z) tables.

Uses published Nbin/area, avoiding an assumed common cosmology. No FKP weights
are applied to counts. This is a published-selection reproduction, not Loa parity.
"""
import argparse
import json
import os
import re
from pathlib import Path
import fitsio
import numpy as np
from mock_selection_test import digest, save

MOCK=Path('/global/cfs/cdirs/desi/public/edr/vac/edr/uchuu/v1.0/BGS-BRIGHT_Uchuu')
DATA=Path('/global/cfs/cdirs/desi/survey/catalogs/edav1/sv3/LSScats')
ALT=Path('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/fuji/LSScats/3.1')
ZE=np.linspace(.15,.55,21)

def table(root, cap):
    p=root/f'BGS_BRIGHT_{cap}_nz.txt'
    d=np.loadtxt(p); sel=(d[:,1]>=.15-1e-8)&(d[:,2]<=.55+1e-8); d=d[sel]
    assert len(d)==20 and np.allclose(d[:,1],ZE[:-1]) and np.allclose(d[:,2],ZE[1:])
    area=float(re.search(r'area is ([0-9.]+)',p.read_text()).group(1))
    return d[:,4]/area,dict(path=str(p),sha256=digest(p),area_deg2=area,
                         statistic='published Nbin per deg2; inherited LSS completeness treatment')

def run(root):
    if (root/'COMPLETE.json').exists(): raise FileExistsError(root)
    h=np.zeros((102,2,20),dtype='i8'); metadata=[]; ranges=[]
    for i in range(102):
        for j,c in enumerate(['S','N']):
            p=MOCK/f'BGS_BRIGHT_{c}_{i}_uchuu.dat.fits'
            stat=[p.stat().st_size,p.stat().st_mtime_ns]
            d=fitsio.read(p,ext=1,columns=['Z','APPMAG_R'])
            good=np.isfinite(d['Z'])&np.isfinite(d['APPMAG_R'])&(d['APPMAG_R']<=19.5)
            h[i,j]=np.histogram(d['Z'][good],bins=ZE)[0]
            ranges.append([float(d['APPMAG_R'].min()),float(d['APPMAG_R'].max())])
            assert stat==[p.stat().st_size,p.stat().st_mtime_ns]
            metadata.append(dict(path=str(p),stat=stat,rows=len(d)))
        if i%20==0: print('Uchuu',i,flush=True)
    original=[]; alternative=[]; tables=[]
    for c in ['S','N']:
        a,m=table(DATA,c);original.append(a);tables.append(m)
        a,m=table(ALT,c);alternative.append(a);tables.append(m)
    original=np.array(original);alternative=np.array(alternative)
    density=h/95.7
    broad=density.reshape(102,2,4,5).sum(-1)
    real=original.reshape(2,4,5).sum(-1)
    np.savez_compressed(root/'HISTOGRAMS.npz',counts=h,z_edges=ZE,
                        data_per_deg2=original,fuji31_per_deg2=alternative)
    save(root/'COMPLETE.json',dict(job=os.environ.get('SLURM_JOB_ID'),
        script_sha256=digest(__file__),histogram_sha256=digest(root/'HISTOGRAMS.npz'),
        mock_sources=metadata,data_sources=tables,magnitude_ranges=ranges,
        caps=['S','N'],shell_edges=[.15,.25,.35,.45,.55],
        mock_to_sv3_mean=(broad.mean(0)/real).tolist(),
        mock_to_sv3_realization_std=(broad.std(0,ddof=1)/real).tolist(),
        fuji31_to_archived=(alternative.reshape(2,4,5).sum(-1)/real).tolist(),
        notes=['Published SV3 selection and area-normalized comparison, not a Loa observed-count comparison.',
               'Raw archived clustering files were not located; inherited Nbin weighting is not independently replayed.',
               'Mock spread is realization standard deviation, not standard error and not independent full-survey covariance.',
               'No FKP weighting, n(z) fitting or catalogue mutation. Geometry and volume are shared between some realizations.',
               'Source for the 95.7 deg2 per hemisphere convention: official Uchuu VAC documentation.']))
    print('COMPLETE Uchuu',broad.mean(0)/real,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(parents=True,exist_ok=True);run(a.root)
