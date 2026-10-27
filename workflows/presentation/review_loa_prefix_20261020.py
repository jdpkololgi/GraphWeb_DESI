"""Reuse the existing 32-draw preflight; no inference or truth access."""
import hashlib
import io
import json
import os
from pathlib import Path
import socket
import tarfile

import numpy as np
from review_loa_preview_20261020 import MODELS, ROOT, THRESHOLD, sha


def main():
    if not socket.gethostname().startswith('nid'):
        raise RuntimeError('Compute node required')
    tile = 'NGC_005_008_003'
    preflight = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_preflight_20261016_v2')
    folder = preflight / 'packed' / (tile + '_32draws')
    meta = json.loads((folder / 'COMPLETE.json').read_text())
    assert sha(folder / 'tile.tar') == meta['archive_sha256']
    assert sha(ROOT / 'INDEX.json') == sha(preflight / 'INDEX.json') == meta['index_sha256']
    result = {}
    with tarfile.open(folder / 'tile.tar') as a, tarfile.open(ROOT / 'packed' / (tile + '_04draws') / 'tile.tar') as b:
        for model in MODELS:
            draws = []
            for d in range(32):
                name = f'{model}/{d:03d}/arrays.npz'
                content = a.extractfile(name).read()
                assert hashlib.sha256(content).hexdigest() == meta['members'][name]['sha256']
                if d < 4:
                    assert content == b.extractfile(name).read()
                with np.load(io.BytesIO(content)) as z:
                    draws.append(z['galaxy_eigen'].astype(float))
            e = np.stack(draws)
            labels = (e > THRESHOLD).sum(axis=-1)
            mean32 = e.mean(axis=0)
            sd32 = e.std(axis=0, ddof=1)
            p32 = np.stack([(labels == k).mean(axis=0) for k in range(4)], axis=-1)
            entries = {}
            for n in (4, 8, 16, 32):
                ep = e[:n]
                pp = np.stack([(labels[:n] == k).mean(axis=0) for k in range(4)], axis=-1)
                delta = ep.mean(axis=0) - mean32
                entries[str(n)] = {
                    'mean_eigen': ep.mean(axis=(0, 1)).tolist(),
                    'tidal_spread': float(((ep[..., 2] - ep[..., 0]) / np.sqrt(2)).mean()),
                    'class_percent': (100 * np.bincount(labels[:n].ravel(), minlength=4) / labels[:n].size).tolist(),
                    'mean_eigen_rmse_to_32': np.sqrt((delta**2).mean(axis=0)).tolist(),
                    'mean_eigen_rms_change_over_sd32': np.sqrt(((delta / np.maximum(sd32, 1e-12))**2).mean(axis=0)).tolist(),
                    'class_probability_tv_to_32': float(.5 * np.abs(pp - p32).sum(axis=-1).mean()),
                }
            groups = [np.bincount(x.ravel(), minlength=4) / x.size for x in np.split(labels, 8)]
            result[model] = {'prefixes': entries,
                             'eight_disjoint_four_draw_class_percent': (100 * np.array(groups)).tolist()}
    output = dict(passed=True, tile=tile, galaxies=len(mean32), first_four_bitwise_equal=True,
                  results=result, archive_sha256=meta['archive_sha256'], source_sha256=sha(__file__),
                  job=os.environ['SLURM_JOB_ID'],
                  caveats=['One pre-existing geometry-selected tile, not survey-wide convergence.',
                           'Prefixes overlap the 32-draw reference; this is not an independent error estimate.',
                           '32 draws are not ground truth; this comparison cannot establish calibration.'])
    with (ROOT / 'PREFIX_DIAGNOSTICS_20261020.json').open('x') as f:
        json.dump(output, f, allow_nan=False)
    print(json.dumps(output), flush=True)


if __name__ == '__main__':
    main()
