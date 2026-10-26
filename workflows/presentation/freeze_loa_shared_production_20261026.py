"""Stage a small immutable, hash-bound source snapshot; never submit jobs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

REPO = Path(__file__).resolve().parents[2]
ILLUSTRIS = REPO.parent / 'TNG/Illustris'
CAMPAIGN = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_shared_sixfit_20261026_v1')
OWN = [
    'configs/loa_shared_production_20261026.json',
    'workflows/presentation/loa_shared_sixfit_20261026.py',
    'workflows/presentation/loa_shared_field_io_20261026.py',
    'workflows/presentation/loa_shared_production_20261026.py',
    'workflows/presentation/loa_shared_production_20261026.slurm',
    'workflows/presentation/test_loa_shared_sixfit_20261026.py',
    'workflows/presentation/test_loa_shared_field_io_20261026.py',
    'workflows/presentation/freeze_loa_shared_production_20261026.py',
]


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 << 20), b''):
            value.update(chunk)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def git(repo, *args):
    return subprocess.check_output(['git', '-C', str(repo), *args])


def freeze():
    root = CAMPAIGN / 'source_production_v1'
    if root.exists():
        raise FileExistsError(root)
    benchmark = read(CAMPAIGN / 'benchmark_ngc_v1/BINDING.json')
    completed = read(CAMPAIGN / 'benchmark_ngc_v1/COMPLETE.json')
    smoke = read(CAMPAIGN / 'smoke_v2/COMPLETE.json')
    if not completed['passed'] or not smoke['passed']:
        raise ValueError('Engineering gates not passed')
    if completed['binding_sha256'] != sha(CAMPAIGN / 'benchmark_ngc_v1/BINDING.json'):
        raise ValueError('Benchmark binding changed')
    sources = {}
    for source, digest in benchmark['source'].items():
        path = Path(source)
        if sha(path) != digest:
            raise ValueError('GPU-tested inference source changed: ' + source)
        base = ILLUSTRIS if path.is_relative_to(ILLUSTRIS) else REPO
        relative = Path('Illustris' if base == ILLUSTRIS else 'GraphWeb_DESI') / path.relative_to(base)
        sources[str(relative)] = path
    for name in OWN:
        committed = git(REPO, 'show', 'HEAD:' + name)
        if hashlib.sha256(committed).hexdigest() != sha(REPO / name):
            raise ValueError('Production source is not committed: ' + name)
        sources['GraphWeb_DESI/' + name] = REPO / name
    # Normalizer coordinate authority and any bounded E2E configuration imports.
    for path in sorted((ILLUSTRIS / 'configs').glob('e2e_*.json')):
        sources['Illustris/' + str(path.relative_to(ILLUSTRIS))] = path
    evidence = {}
    for name in ('smoke_v2/COMPLETE.json', 'smoke_v2/BINDING.json',
                 'benchmark_ngc_v1/COMPLETE.json', 'benchmark_ngc_v1/BINDING.json'):
        path = CAMPAIGN / name
        evidence[name] = dict(path=str(path), sha256=sha(path))
    root.mkdir()
    records = {}
    provenance = {}
    for relative, path in sorted(sources.items()):
        destination = root / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        digest = sha(path)
        shutil.copyfile(path, destination)
        if sha(destination) != digest or sha(path) != digest:
            raise ValueError('Source changed while staging: ' + str(path))
        destination.chmod(0o444)
        records[relative] = digest
        provenance[relative] = str(path)
    value = dict(schema='loa-shared-source-freeze-v1', snapshot_root=str(root),
        snapshot_files=records, original_sources=provenance, evidence=evidence,
        input_index_sha256=benchmark['index_sha256'],
        git={name: git(repo, 'rev-parse', 'HEAD').decode().strip()
             for name, repo in (('GraphWeb_DESI', REPO), ('Illustris', ILLUSTRIS))},
        inference_sources_identical_to_gpu_benchmark=True,
        new_layer='density writer, independent storage/coverage checks and batch launch only',
        automatic_retries=False, posterior_qualified=False)
    with (root / 'FREEZE.json').open('x') as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
    (root / 'FREEZE.json').chmod(0o444)
    (CAMPAIGN / 'production_logs_v1').mkdir(exist_ok=False)
    print(json.dumps(dict(snapshot=str(root), manifest_sha256=sha(root / 'FREEZE.json'),
                          files=len(records), git=value['git'])))


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    freeze()
