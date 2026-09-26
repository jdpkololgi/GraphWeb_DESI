"""Read-only Loa source revalidation; does not select galaxies or run inference."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def audit(registry):
    if not os.environ.get('SLURM_JOB_ID') or not socket.gethostname().startswith('nid'):
        raise RuntimeError('live catalogue hashing requires a compute allocation')
    import fitsio
    import numpy as np
    registry_sha = digest(registry)
    sources = json.loads(registry.read_text())['desi_candidate']
    entries = [('data_'+k, v) for k,v in sources['data'].items()]
    entries += [(kind+'_'+str(i),v) for kind in ('full_random','clustering_random')
                for i,v in enumerate(sources[kind])]
    records = []
    for kind, expected in entries:
        path = Path(expected['path'])
        before = path.stat()
        with fitsio.FITS(path) as f:
            rows = f[1].get_nrows()
            columns = f[1].get_colnames()
            checks = dict(size=before.st_size == expected['bytes'],
                          rows=rows == expected['rows'], columns=columns == expected['columns'])
            sample = None
            if kind == 'data_full':
                wanted = ['TARGETID','RA','DEC','Z_not4clus','ZWARN','DELTACHI2',
                          'SPECTYPE','GOODHARDLOC','GOODPRI','LOCATION_ASSIGNED',
                          'FRAC_TLOBS_TILES','FRACZ_TILELOCID','mod_success_rate']
                index = np.unique(np.linspace(0, rows-1, 8192, dtype=np.int64))
                data = f[1].read(rows=index, columns=wanted)
                sample = dict(rows=index.tolist(), columns=wanted,
                              sha256=hashlib.sha256(data.tobytes()).hexdigest(),
                              finite_counts={k:int(np.isfinite(data[k]).sum()) for k in wanted
                                             if data[k].dtype.kind in 'if'},
                              scope='schema/finite-value diagnostic; no success selection adopted')
        live_hash = digest(path) if kind.startswith('data_') else None
        if live_hash is not None:
            checks['sha256'] = live_hash == expected['sha256']
        after = path.stat()
        checks['stable_during_read'] = (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
        records.append(dict(kind=kind,path=str(path),bytes=before.st_size,mtime_ns=before.st_mtime_ns,
                            rows=rows,columns=columns,checks=checks,sha256=live_hash,
                            content_hash_verified=live_hash is not None,sample=sample))
    if digest(registry) != registry_sha:
        raise RuntimeError('registry changed during audit')
    return dict(schema='p12a-loa-source-audit-v1',registry=str(registry),registry_sha256=registry_sha,
                source_sha256=digest(__file__),job=os.environ['SLURM_JOB_ID'],node=socket.gethostname(),
                records=records,pass_checks=all(all(r['checks'].values()) for r in records),
                random_content_hashes_verified=False,success_policy_frozen=False,
                ready_for_desi_canary=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--registry',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    result = audit(args.registry)
    with args.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False); f.write('\n')
    print(json.dumps(dict(output=str(args.output),pass_checks=result['pass_checks'])))
    raise SystemExit(0 if result['pass_checks'] else 1)
