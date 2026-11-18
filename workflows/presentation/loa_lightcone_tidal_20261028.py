"""Full-cap free-space T-web display layers from existing R7 density draws.

CPU allocation only. Unknown exterior (including the opposite cap) has delta=0.
Classifications are exact on the declared threshold grid, before projection.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import socket
import sys
import time

import h5py
import numpy as np
from scipy import fft

sys.path.insert(0, '/global/u2/d/dkololgi/TNG/Illustris/workflows/sbi')
import e2e_footprint_tidal_operator as op
from loa_lightcone_20261028 import BASE, north, sha


def check_operator():
    rng = np.random.default_rng(17)
    delta = rng.normal(size=(4, 3, 5))
    k = op.kernel(delta.shape)
    fast = op.convolve(delta, k, workers=4)
    direct = np.zeros_like(fast)
    for t in np.ndindex(delta.shape):
        for s in np.ndindex(delta.shape):
            direct[t] += delta[s] * op.cube_tensor(*[t[i]-s[i] for i in range(3)])
    errors = dict(direct=float(abs(fast-direct).max()),
                  padding=float(abs(fast-op.convolve(delta, op.kernel(delta.shape, extra_padding=3), workers=4)).max()),
                  trace=float(abs(fast[..., 0]+fast[..., 3]+fast[..., 5]-delta).max()),
                  self=float(abs(op.cube_tensor(0, 0, 0)-np.array([1/3,0,0,1/3,0,1/3])).max()))
    assert max(errors.values()) < 1e-10, errors
    return errors


def main(args):
    assert os.environ.get('SLURM_JOB_ID') and socket.gethostname().startswith('nid')
    out = args.output
    out.mkdir(parents=True, exist_ok=False)
    manifest = json.loads(args.manifest.read_text())
    thresholds = np.round(np.arange(-.2, 1.0001, .05), 2)
    assets = {}
    record = dict(schema='loa-lightcone-tidal-v1', thresholds=thresholds.tolist(),
        encoding='uint8 [draw, radial, RA, eigenvalue]; code = count(thresholds < eigenvalue); 255 = missing',
        draws=4, density_manifest_sha256=sha(args.manifest), caps={}, assets=assets,
        operator='Cuboid-integrated free-space potential Hessian; Laplacian(Phi)=delta; no added smoothing',
        source='Within each cap independently: delta=(rho-1)*fractional random support where prediction coverage>0; delta=0 elsewhere, including the opposite cap and all unmodelled exterior. No periodic images or box mean subtraction.',
        limitations='Exploratory closure-dependent T-web, not calibrated real-DESI environments or a connectivity measurement. Fractional footprint occupancy is a source assumption, not completeness correction. Projected class fractions combine angular samples, and optionally four draws; they are not calibrated probabilities. Mixed/tied dominant classes are shown separately.',
        provenance=dict(script_sha256=sha(__file__), operator_sha256=sha(op.__file__),
                        job=os.environ['SLURM_JOB_ID'], node=socket.gethostname()))
    qa = dict(passed=False, operator=check_operator(), caps=[])
    re = np.array(manifest['radial_edges']); rc = (re[1:]+re[:-1])/2
    ac = (np.array(manifest['ra_edges'])[1:]+np.array(manifest['ra_edges'])[:-1])/2
    n = len(rc)*len(ac)
    start = time.monotonic()
    for cap, capvalue in [('NGC', True), ('SGC', False)]:
        source = BASE/f'{cap}_seed17/unet_x_17_FIELD.h5'
        assert sha(source) == manifest['caps'][cap]['field_sha256']
        slices = manifest['caps'][cap]['slices']
        with h5py.File(source) as f, h5py.File(out/f'{cap}-sampled-eigenvalues.h5','w') as stored:
            shape = f['coverage'].shape; origin = np.array(f.attrs['origin_mpc_h'])
            coverage = f['coverage'][:]; support = f['support'][:]
            assert f['rho'].shape == (4,*shape) and f.attrs['cell_mpc_h'] == 6.766
            locations=[]; validities=[]
            for item in slices:
                d = np.deg2rad(np.mean(item['dec'])); a=np.deg2rad(ac)[None,:]; r=rc[:,None]
                xyz=np.stack(np.broadcast_arrays(r*np.cos(d)*np.cos(a),r*np.cos(d)*np.sin(a),r*np.sin(d)),axis=-1)
                idx=np.floor((xyz-origin)/6.766).astype(int)
                valid=np.all((idx>=0)&(idx<np.array(shape)),axis=-1)&(north(ac,np.mean(item['dec']))[None,:]==capvalue)
                clipped=np.clip(idx,0,np.array(shape)-1)
                loc=tuple(clipped[...,j] for j in range(3))
                valid &= (coverage[loc]>0)&(support[loc]>0)
                # Verify agreement with the already-published density support.
                old=np.frombuffer(gzip.decompress((args.manifest.parent/item['density']).read_bytes()),dtype='<f4').reshape(4,len(rc),len(ac))
                assert np.array_equal(valid,np.all(old>0,axis=0))
                locations.append(np.ravel_multi_index(loc,shape).reshape(-1))
                validities.append(valid.reshape(-1))
                stored.create_dataset(str(item['bin']), (4,n,3), dtype='<f4',compression='gzip')
            print(json.dumps(dict(cap=cap,event='kernel_start',shape=shape)),flush=True)
            kernels=op.kernel(shape)
            codes=[np.full((4,n,3),255,dtype='u1') for _ in slices]
            trace_max=0.; exact=True
            for draw in range(4):
                delta=op.source_delta(f['rho'][draw],support,coverage)
                spectrum=fft.rfftn(delta,s=kernels.shape[1:],workers=args.workers)
                tensor=np.empty((*shape,6),np.float64)
                crop=tuple(slice(0,x) for x in shape)
                for c in range(6):
                    ks=fft.rfftn(kernels[c],workers=args.workers)
                    ks *= spectrum
                    tensor[...,c]=fft.irfftn(ks,s=kernels.shape[1:],workers=args.workers)[crop]
                    del ks
                del spectrum
                err=float(abs(tensor[...,0]+tensor[...,3]+tensor[...,5]-delta).max())
                assert err<1e-8, err
                trace_max=max(trace_max,err)
                flat=tensor.reshape(-1,6)
                for j,item in enumerate(slices):
                    valid=validities[j]
                    ev=np.full((n,3),np.nan)
                    ev[valid]=op.eigensystem(flat[locations[j][valid]])
                    assert np.isfinite(ev[valid]).all() and np.all(np.diff(ev[valid],axis=-1)>=0)
                    # searchsorted(left) preserves strict > even on threshold ties.
                    code=np.searchsorted(thresholds,ev[valid],side='left').astype('u1')
                    codes[j][draw,valid]=code
                    stored[str(item['bin'])][draw]=ev.astype('<f4')
                    for q,threshold in enumerate(thresholds):
                        assert np.array_equal(np.sum(code>q,axis=-1),np.sum(ev[valid]>threshold,axis=-1))
                del tensor,delta
                print(json.dumps(dict(cap=cap,draw=draw,trace_error=err,seconds=time.monotonic()-start)),flush=True)
            del kernels
            entries=[]
            for j,item in enumerate(slices):
                name=f'{cap}-{item["bin"]}-tidal.bin'
                raw=codes[j].tobytes(); packed=gzip.compress(raw,compresslevel=9,mtime=0)
                assert gzip.decompress(packed)==raw
                (out/name).write_bytes(packed)
                assets[name]=dict(bytes=len(packed),sha256=hashlib.sha256(packed).hexdigest())
                entry=dict(bin=item['bin'],file=name,valid_pixels=int(validities[j].sum()))
                if validities[j].any():
                    pixel=int(np.flatnonzero(validities[j])[len(np.flatnonzero(validities[j]))//2])
                    entry['reference']=dict(pixel=pixel,codes=codes[j][:,pixel,:].tolist(),eigenvalues=stored[str(item['bin'])][:,pixel,:].tolist())
                entries.append(entry)
            record['caps'][cap]=dict(slices=entries,field_sha256=manifest['caps'][cap]['field_sha256'])
            qa['caps'].append(dict(cap=cap,trace_max=trace_max,threshold_classification_exact=exact,support_matches_density=True))
        del codes,coverage,support,locations,validities
    qa.update(passed=True,asset_bytes=sum(x['bytes'] for x in assets.values()),job=os.environ['SLURM_JOB_ID'],seconds=time.monotonic()-start)
    assert qa['asset_bytes']<38*1024**2, 'Tidal payload exceeds remaining site budget; preserve outputs and review packaging'
    (out/'manifest.json').write_text(json.dumps(record,separators=(',',':')))
    (out/'QA.json').write_text(json.dumps(qa,indent=2))
    print(json.dumps(qa),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--workers',type=int,default=32)
    main(p.parse_args())
