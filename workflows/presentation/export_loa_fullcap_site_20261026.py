"""Export all saved Loa tile volumes for the private six-fit 3D website.

No inference; no cross-tile interpolation; no individual identifiers exported.
Display values pool 2x2x2 native voxels, matching the older viewer resolution.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import socket
import tarfile
import time
import numpy as np

ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1')
MODELS = [f'{f}_{s}' for f in ('legacy_v', 'unet_x', 'transformer_v') for s in (17, 29)]
NAMES = ['Earlier CFM · seed 17', 'Earlier CFM · seed 29', 'Expanded U-Net · seed 17', 'Expanded U-Net · seed 29', 'Transformer · seed 17', 'Transformer · seed 29']

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()

def read(path):
    return json.loads(path.read_text())

def pool(a):
    assert a.shape == (32, 16, 16)
    return a.reshape(16, 2, 8, 2, 8, 2).mean(axis=(1, 3, 5))

def quantize(a, supported):
    assert np.isfinite(a).all() and (a > 0).all()
    logs = np.log2(a)
    assert logs.min() >= -8 and logs.max() <= 8
    q = np.rint((logs + 8) * (65534 / 16)).astype('<u2')
    error = float(np.max(np.abs(2 ** (q.astype(float) * (16 / 65534) - 8) / a - 1)))
    assert error < .000085
    q[:, ~supported] = 65535
    return q, error

def extract(case, audit, sample):
    directory = ROOT/'packed'/(case['id']+'_04draws')
    meta = read(directory/'COMPLETE.json')
    assert meta['passed'] and sha(directory/'tile.tar') == audit['archive_sha256'] == meta['archive_sha256']
    with tarfile.open(directory/'tile.tar') as tar:
        def load(member, columns):
            raw = tar.extractfile(member).read()
            assert hashlib.sha256(raw).hexdigest() == meta['members'][member]['sha256']
            with np.load(io.BytesIO(raw)) as f:
                return {k: f[k] for k in columns}
        prep = load('prepared/arrays.npz', ['old_joint', 'xyz', 'Z'])
        assert len(prep['xyz']) == audit['retained']
        support = pool(prep['old_joint'][1, 16:48, 16:32, 16:32].astype(float)) >= .5
        density = []
        error = 0.
        for model in MODELS:
            draws = np.stack([pool(load(f'{model}/{d:03d}/arrays.npz', ['rho'])['rho'].astype(float)) for d in range(4)])
            values = np.concatenate([draws.mean(axis=0, keepdims=True), draws])
            q, e = quantize(values, support)
            density.append(q)
            error = max(error, e)
        points = np.column_stack([prep['xyz'][sample], prep['Z'][sample]]).astype('<f4')
        tile = {'id':case['id'], 'origin':case['core_origin_mpc_h'], 'key':case['key'],
                'median_z':case['z_median'], 'shell':case['shell'], 'n':audit['retained'], 'eligible':audit['eligible'],
                'counts':[audit['models'][m]['retained_class_draw_counts'] for m in MODELS]}
        return tile, np.stack(density), points, error, int(support.sum())

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert socket.gethostname().startswith('nid') and os.environ.get('SLURM_JOB_ID')
    start = time.monotonic()
    index, audit = read(ROOT/'INDEX.json'), read(ROOT/'PREVIEW_AUDIT.json')
    assert audit['passed'] and audit['tiles'] == len(index['cases']) == 2472
    by = {r['tile']:r for r in audit['per_tile']}
    args.output.mkdir(exist_ok=False)
    manifest = {'schema':1, 'kind':'mosaic', 'models':MODELS, 'names':NAMES, 'draws':4,
                'shape':[16,8,8], 'spacing':13.532, 'native_spacing':6.766,
                'quantization':{'log2_min':-8,'step':16/65534,'missing':65535}, 'caps':{}}
    receipt = {'passed':False, 'job':os.environ['SLURM_JOB_ID'], 'node':socket.gethostname(),
               'source_sha256':sha(Path(__file__)), 'index_sha256':sha(ROOT/'INDEX.json'),
               'audit_sha256':sha(ROOT/'PREVIEW_AUDIT.json'), 'assets':{}, 'caps':{}}
    def save(name, raw):
        path = args.output/name
        compressed = gzip.compress(raw, compresslevel=6, mtime=0)
        path.write_bytes(compressed)
        assert gzip.decompress(path.read_bytes()) == raw
        receipt['assets'][name] = {'bytes':len(compressed), 'sha256':sha(path)}
        return 'six-fit-data/'+name
    for cap in ('NGC','SGC'):
        cases = [c for c in index['cases'] if c['cap']==cap]
        total = sum(by[c['id']]['retained'] for c in cases)
        selected = np.linspace(0, total-1, min(total,100000), dtype=np.int64)
        offset = 0
        work = []
        for c in cases:
            n = by[c['id']]['retained']
            sl = selected[(selected>=offset)&(selected<offset+n)]-offset
            work.append((c,by[c['id']],sl))
            offset += n
        results = []
        with ThreadPoolExecutor(max_workers=8) as executor:
            for i, r in enumerate(executor.map(lambda x:extract(*x), work),1):
                results.append(r)
                if i%200==0:
                    print(cap, i, 'tiles', round(time.monotonic()-start,1),'s',flush=True)
        assert len(results)==len(cases)
        data = np.stack([r[1] for r in results])
        assert data.shape == (len(cases),6,5,16,8,8)
        tiles = [r[0] for r in results]
        points = np.concatenate([r[2] for r in results])
        assert len(points)==min(total,100000)
        origins = np.array([t['origin'] for t in tiles])
        bounds = {k:[float(origins[:,j].min()),float((origins[:,j]+[32,16,16][j]*6.766).max())] for j,k in enumerate('xyz')}
        files = [[save(f'{cap}-{model}-{d}.bin',data[:,m,d].astype('<u2').tobytes()) for d in range(5)] for m,model in enumerate(MODELS)]
        manifest['caps'][cap] = {'tiles':tiles, 'bounds':bounds, 'retained':total,
                                 'eligible':sum(t['eligible'] for t in tiles),
                                 'files':files,'galaxies':save(cap+'-galaxies.bin',points.tobytes()),
                                 'galaxy_count':len(points)}
        receipt['caps'][cap] = {'tiles':len(tiles),'retained':total,'context_galaxies':len(points),
                                'supported_cells':sum(r[4] for r in results),
                                'max_relative_density_error':max(r[3] for r in results)}
        print(cap, receipt['caps'][cap], flush=True)
    manifest['note'] = 'Independent local-tile posteriors, not joint full-cap draws. 2× mean pooling in all three dimensions; no point thinning or cross-tile interpolation. Pooled support >=0.5. Mean and four individual draws for six frozen fits. Context galaxies are a deterministic sample of retained positions, not environment-labelled galaxies.'
    (args.output/'manifest.json').write_text(json.dumps(manifest,separators=(',',':'),allow_nan=False))
    receipt.update(passed=True, seconds=time.monotonic()-start,manifest_sha256=sha(args.output/'manifest.json'),
                   total_bytes=sum(x['bytes'] for x in receipt['assets'].values()),global_joint_posterior=False)
    (args.output/'PACKAGING.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps({k:v for k,v in receipt.items() if k!='assets'}),flush=True)

if __name__ == '__main__':
    main()
