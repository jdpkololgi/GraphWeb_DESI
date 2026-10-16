"""Recover Gaussian smoothing from saved TNG density and tidal trace.

One sequential read of 4 GiB in small slabs. Plane sums reduce the 3D Gaussian
transfer to three 1D Fourier tests; no new 3D FFT or model inference.
"""
import json
import struct
import zipfile
from pathlib import Path

import numpy as np

OUT = Path('/global/homes/d/dkololgi/fmsu_tokyo_20261010/legacy_selection_tests_20261014')
DENS = Path('/pscratch/sd/d/dkololgi/IllustrisTNG300_densities/snap_099_total_dens_0.npz')
EIG = Path('/pscratch/sd/d/dkololgi/Cosmic_env_TNG_cache/snap_099_eig_vals_0.npz')
CWEB = Path('/pscratch/sd/d/dkololgi/Cosmic_env_TNG_cache/TNG300_snap_099_tweb_env_merged.npz')


def member(path, name):
    with zipfile.ZipFile(path) as z:
        info=z.getinfo(name+'.npy'); assert info.compress_type == zipfile.ZIP_STORED
    with path.open('rb') as f:
        f.seek(info.header_offset); header=f.read(30)
        assert header[:4] == b'PK\x03\x04'
        fn,extra=struct.unpack_from('<HH',header,26);f.seek(fn+extra,1)
        version=np.lib.format.read_magic(f)
        reader=np.lib.format.read_array_header_1_0 if version==(1,0) else np.lib.format.read_array_header_2_0
        shape,fortran,dtype=reader(f); assert not fortran
        return dict(path=path,shape=shape,dtype=dtype,offset=f.tell())


def slab(m, start, count):
    with m['path'].open('rb') as f:
        f.seek(m['offset']+start*m['dtype'].itemsize)
        return np.fromfile(f,dtype=m['dtype'],count=count)


def main():
    OUT.mkdir(exist_ok=True)
    d=member(DENS,'dens');e=member(EIG,'eig_vals');c=member(CWEB,'cweb');n=d['shape'][0]
    assert d['shape']==(n,n,n) and e['shape']==(3,n,n,n) and c['shape']==(n,n,n)
    with np.load(DENS) as f: box=float(f['Boxsize'])
    density=np.zeros((3,n));trace=np.zeros((3,n))
    for start in range(0,n,8):
        width=min(8,n-start); dd=slab(d,start*n*n,width*n*n).reshape(width,n,n)
        tt=np.zeros_like(dd)
        for j in range(3):tt += slab(e,j*n**3+start*n*n,width*n*n).reshape(dd.shape)
        for a,projection in [(dd,density),(tt,trace)]:
            projection[0,start:start+width]=a.sum((1,2))
            projection[1]+=a.sum((0,2));projection[2]+=a.sum((0,1))
    mean=density[0].sum()/n**3
    result={'boxsize':box,'ngrid':n,'files':[str(DENS),str(EIG),str(CWEB)],'axes':[]}
    for axis in range(3):
        dk=np.fft.rfft(density[axis]/n**2/mean-1);tk=np.fft.rfft(trace[axis]/n**2)
        modes=np.arange(1,26);k=2*np.pi*modes/box;ratio=tk[modes]/dk[modes]
        slope,intercept=np.polyfit(k*k,np.log(np.abs(ratio)),1)
        fitted=np.exp(intercept+slope*k*k)
        result['axes'].append(dict(axis=axis,gaussian_radius=float(np.sqrt(max(0.,-2*slope))),
                                  transfer_amplitude=float(np.exp(intercept)),
                                  max_relative_fit_error=float(np.max(np.abs(np.abs(ratio)/fitted-1))),
                                  max_transfer_phase_radians=float(np.max(np.abs(np.angle(ratio)))),
                                  modes=modes.tolist()))
    rng=np.random.default_rng(20261014);flat=rng.choice(n**3,8192,replace=False)
    mm=np.memmap(EIG,mode='r',dtype=e['dtype'],offset=e['offset'],shape=(3,n**3))
    eig=np.asarray(mm[:,flat]);cm=np.memmap(CWEB,mode='r',dtype=c['dtype'],offset=c['offset'],shape=(n**3,))
    labels=np.asarray(cm[flat]);unique=np.unique(labels)
    result['class_values']=unique.tolist()
    result['threshold_checks']={str(t):int(np.count_nonzero((eig>t).sum(0)!=labels)) for t in [0.,.1,.2,.25,.3,.4,.5,.6,1.]}
    result['sampled_voxels']=len(flat)
    result['interpretation']='Radius uses the saved box-length units. Exact Gaussian transfer plus threshold parity binds these saved fields; it does not alone prove which labels trained a checkpoint.'
    (OUT/'label_scale.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
