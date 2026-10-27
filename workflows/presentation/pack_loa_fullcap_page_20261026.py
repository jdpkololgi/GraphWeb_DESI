"""Package private, reduced-resolution independent-tile mosaics; no inference."""
import base64
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
OUT = ROOT / 'page_fullcap_20261026_v2'
MODELS = [f'{f}_{s}' for f in ('legacy_v', 'unet_x', 'transformer_v') for s in (17, 29)]
NAMES = ['Earlier CFM 17', 'Earlier CFM 29', 'Expanded U-Net 17', 'Expanded U-Net 29', 'Transformer 17', 'Transformer 29']

def read(p):
    return json.loads(p.read_text())

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()

def encode(v):
    return base64.b64encode(gzip.compress(v, compresslevel=9, mtime=0)).decode()

def tile(case, audit):
    folder = ROOT / 'packed' / (case['id'] + '_04draws')
    meta = read(folder / 'COMPLETE.json')
    assert meta['passed'] and meta['archive_sha256'] == audit['archive_sha256']
    assert sha(folder / 'tile.tar') == audit['archive_sha256']
    with tarfile.open(folder / 'tile.tar') as t:
        def load(name, keys):
            raw = t.extractfile(name).read()
            assert hashlib.sha256(raw).hexdigest() == meta['members'][name]['sha256']
            with np.load(io.BytesIO(raw)) as f:
                return {k: f[k] for k in keys}
        joint = load('prepared/arrays.npz', ['old_joint'])['old_joint']
        # Legacy native core starts at [16,16,16]; central z plane is index 8.
        support = joint[1, 16:48, 16:32, 24]
        assert support.shape == (32, 16)
        mask = (support.reshape(8, 4, 4, 4).min(axis=(1, 3)) >= .5)
        values = []
        for m in MODELS:
            planes = []
            for d in range(4):
                rho = load(f'{m}/{d:03d}/arrays.npz', ['rho'])['rho']
                assert rho.shape == (32, 16, 16) and np.isfinite(rho).all() and (rho > 0).all()
                planes.append(rho[:, :, 8].astype(float).reshape(8, 4, 4, 4).mean(axis=(1, 3)))
            values.append(np.mean(planes, axis=0))
        v = np.array(values)
        # 16 logarithmic display levels, common scale for all six fits and caps.
        logs = np.log2(v)
        assert logs.min() >= -4 and logs.max() <= 4
        q = np.rint((logs + 4) * 15 / 8).astype('u1').ravel()
        packed = (q[::2] | (q[1::2] << 4)).tobytes()
        np.testing.assert_array_equal(np.stack([np.frombuffer(packed, 'u1') & 15, np.frombuffer(packed, 'u1') >> 4], axis=1).ravel(), q)
        error = float(np.max(np.abs(2 ** (q.reshape(v.shape) * (8 / 15) - 4) / v - 1)))
        return dict(id=case['id'], key=case['key'], origin=case['core_origin_mpc_h'], shell=case['shell'],
                    z=case['z_median'], n=audit['retained'], eligible=audit['eligible'],
                    mask=np.packbits(mask.ravel(), bitorder='little').tolist(),
                    counts=[audit['models'][m]['retained_class_draw_counts'] for m in MODELS],
                    density=packed, error=error, valid=int(mask.sum()))

def main():
    assert socket.gethostname().startswith('nid') and os.environ.get('SLURM_JOB_ID')
    start = time.monotonic()
    index, audit = read(ROOT / 'INDEX.json'), read(ROOT / 'PREVIEW_AUDIT.json')
    assert audit['passed'] and audit['tiles'] == 2472
    review = read(ROOT / 'SCIENTIFIC_DIAGNOSTICS_20261020.json')
    assert review['passed'] and review['audit_sha256'] == sha(ROOT / 'PREVIEW_AUDIT.json')
    OUT.mkdir(exist_ok=False)
    by = {x['tile']: x for x in audit['per_tile']}
    records = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        for i, r in enumerate(pool.map(lambda c: tile(c, by[c['id']]), index['cases']), 1):
            records.append(r)
            if i % 200 == 0:
                print(i, 'tiles', round(time.monotonic()-start, 1), 'seconds', flush=True)
    assert len(records) == len({r['id'] for r in records}) == 2472
    assert sum(r['n'] for r in records) == audit['retained']
    template = Path(__file__).with_name('loa_fullcap_page_20261026.html').read_text()
    receipt = {'passed': True, 'job': os.environ['SLURM_JOB_ID'], 'node': socket.gethostname(),
               'sources': {n: sha(ROOT/n) for n in ('INDEX.json', 'PREVIEW_AUDIT.json', 'SCIENTIFIC_DIAGNOSTICS_20261020.json', 'PREFIX_DIAGNOSTICS_20261020.json')},
               'code_sha256': sha(Path(__file__)), 'template_sha256': sha(Path(__file__).with_name('loa_fullcap_page_20261026.html')),
               'display': 'Four-draw mean; native central z-plane 8; xy 4x4 arithmetic pooling; conservative min support >=0.5; 16 log2 density levels over [-4,4]; no interpolation between tiles',
               'caps': {}}
    for cap in ('NGC', 'SGC'):
        rows = [r for r in records if r['id'].startswith(cap)]
        # Private source IDs and individual galaxy properties are never exported.
        data = {'cap': cap, 'models': MODELS, 'names': NAMES, 'cell': 6.766, 'tileCount': len(rows),
                'retained': sum(r['n'] for r in rows), 'eligible': sum(r['eligible'] for r in rows),
                'tiles': [[r['key'], [round(x, 2) for x in r['origin']], round(r['z'], 4), r['n'], r['eligible'], r['mask'], r['counts'], r['shell']] for r in rows],
                'density': base64.b64encode(b''.join(r['density'] for r in rows)).decode(),
                'summary': review['summary'][cap+'_all']['models']}
        for j, m in enumerate(MODELS):
            np.testing.assert_array_equal(np.sum([r['counts'][j] for r in rows], axis=0), data['summary'][m]['class_counts'])
        packed = encode(json.dumps(data, separators=(',', ':'), allow_nan=False).encode())
        html = template.replace('__PAYLOAD__', packed)
        assert len(html.encode()) <= 256*1024, (cap, len(html.encode()))
        path = OUT / (cap + '.html')
        path.write_text(html)
        receipt['caps'][cap] = {'tiles': len(rows), 'retained': data['retained'], 'eligible': data['eligible'],
                               'empty_retained_tiles': sum(r['n']==0 for r in rows),
                               'supported_display_cells': sum(r['valid'] for r in rows),
                               'max_relative_quantization_error': max(r['error'] for r in rows),
                               'bytes': path.stat().st_size, 'sha256': sha(path)}
    receipt['seconds'] = time.monotonic()-start
    (OUT / 'PACKAGING.json').write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt), flush=True)

if __name__ == '__main__':
    main()
