"""Bounded production/archive lifecycle; never submit, retry or requeue jobs.

Only redundant, checksum-matched draw files under this campaign's literal
ROOT may be removed. Prepared inputs, archives and preflight outputs survive.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from types import SimpleNamespace

import loa_survey_runner_20261016 as r
import pack_loa_survey_tile_20261016 as packing

ROOT = Path('/pscratch/sd/d/dkololgi/graphweb_desi/outputs/loa_survey_production_20261016_v1')
PREFLIGHT = r.DEFAULT_ROOT
HERE = Path(__file__).resolve().parent
PYTHON = '/pscratch/sd/d/dkololgi/conda/envs/cosmic_env/bin/python'
CAMPAIGN = 'loa-six-fit-20261016-v1'


def expected_members(draws):
    if draws not in (4, 32):
        raise ValueError('Only registered ensemble sizes')
    bundles = ['prepared']
    for seed in r.SEEDS:
        for draw in range(draws):
            bundles.append(f'coarse_{seed}/{draw:03d}')
            bundles.extend(f'{a}_{p}_{seed}/{draw:03d}' for a, p in r.FAMILIES)
    return {f'{bundle}/{name}' for bundle in bundles for name in ('arrays.npz', 'COMPLETE.json')}


def archive_path(root, tile, draws):
    if not tile or PurePosixPath(tile).name != tile or tile in ('.', '..'):
        raise ValueError('Unsafe tile identifier')
    return root/'packed'/f'{tile}_{draws:02d}draws'


def verify_archive(root, tile, draws):
    """Verify container, exact whitelist, every payload and scientific bindings."""
    folder = archive_path(root, tile, draws)
    meta = r.read(folder/'COMPLETE.json')
    index = r.read(root/'INDEX.json')
    if (not meta['passed'] or meta['tile'] != tile or meta['draws_per_fit'] != draws
            or meta['index_sha256'] != r.sha(root/'INDEX.json')
            or tile not in {x['id'] for x in index['cases']}
            or set(meta['members']) != expected_members(draws)
            or r.sha(folder/'tile.tar') != meta['archive_sha256']):
        raise ValueError('Archive identity/hash/membership mismatch')
    receipts = {}
    with tarfile.open(folder/'tile.tar', 'r') as archive:
        members = archive.getmembers()
        if len(members) != len(meta['members']) or {x.name for x in members} != set(meta['members']):
            raise ValueError('Duplicate/missing/unexpected tar members')
        for member in members:
            if not member.isfile() or member.size != meta['members'][member.name]['bytes']:
                raise ValueError('Unsafe tar entry/type/size')
            digest = hashlib.sha256()
            with archive.extractfile(member) as stream:
                if member.name.endswith('COMPLETE.json'):
                    data = stream.read()
                    receipts[member.name] = json.loads(data)
                    digest.update(data)
                else:
                    for block in iter(lambda: stream.read(8 << 20), b''):
                        digest.update(block)
            if digest.hexdigest() != meta['members'][member.name]['sha256']:
                raise ValueError('Archive member hash mismatch')
    for name, receipt in receipts.items():
        payload = name.removesuffix('COMPLETE.json')+'arrays.npz'
        if not receipt['passed'] or receipt['sha256'] != meta['members'][payload]['sha256']:
            raise ValueError('Receipt/payload mismatch')
        if receipt.get('tile') != tile:
            raise ValueError('Wrong tile in receipt')
    prepared = receipts['prepared/COMPLETE.json']
    if prepared['binding'] != r.digest(index):
        raise ValueError('Prepared/index binding mismatch')
    for seed in r.SEEDS:
        for draw in range(draws):
            coarse = receipts[f'coarse_{seed}/{draw:03d}/COMPLETE.json']
            if (coarse['seed'], coarse['draw']) != (seed, draw):
                raise ValueError('Coarse address mismatch')
            if coarse['checkpoint_sha256'] != index['checkpoints'][f'coarse_{seed}']['sha256']:
                raise ValueError('Coarse checkpoint mismatch')
            for architecture, prediction in r.FAMILIES:
                model = f'{architecture}_{prediction}_{seed}'
                fine = receipts[f'{model}/{draw:03d}/COMPLETE.json']
                ck = f'fine_{seed}' if architecture == 'legacy' else model
                if ((fine['seed'], fine['draw'], fine['model']) != (seed, draw, model)
                        or fine['binding'] != r.digest([coarse['binding'], model, draw, coarse['sha256']])
                        or fine['coarse_sha256'] != coarse['sha256']
                        or fine['checkpoint_sha256'] != index['checkpoints'][ck]['sha256']
                        or fine['truth_opened'] is not False):
                    raise ValueError('Fine address/checkpoint/coarse binding mismatch')
    return meta


def safe_destination(base, name):
    parts = PurePosixPath(name).parts
    if PurePosixPath(name).is_absolute() or not parts or any(p in ('.', '..') for p in parts):
        raise ValueError('Unsafe archive path')
    destination = base.joinpath(*parts)
    for candidate in (base, *destination.parents, destination):
        if candidate.is_symlink():
            raise ValueError('Symlink in staging path')
    if not destination.resolve().is_relative_to(base.resolve()):
        raise ValueError('Escaping staging path')
    return destination


def restore(root, tile, draws):
    """Fill missing files, verify existing ones; never overwrite or extract links."""
    meta = verify_archive(root, tile, draws)
    base = root/'tiles'/tile
    with tarfile.open(archive_path(root, tile, draws)/'tile.tar', 'r') as archive:
        for member in archive:
            destination = safe_destination(base, member.name)
            expected = meta['members'][member.name]['sha256']
            if destination.exists():
                if not destination.is_file() or r.sha(destination) != expected:
                    raise ValueError('Existing staged file differs; preserved')
                continue
            destination.parent.mkdir(parents=True, exist_ok=True)
            fd, temporary = tempfile.mkstemp(prefix='.restore.', dir=destination.parent)
            with os.fdopen(fd, 'wb') as output, archive.extractfile(member) as source:
                shutil.copyfileobj(source, output)
                output.flush(); os.fsync(output.fileno())
            if r.sha(temporary) != expected:
                raise ValueError('Restored payload mismatch')
            os.link(temporary, destination)
            os.unlink(temporary)
    return meta


def cleanup(root, tile, draws):
    """Remove only whitelisted redundant draw files, after full archive checking.

    Not a general recursive cleanup. An interruption can be resumed safely;
    prepared files, unknown files and hidden partial bundles are never deleted.
    """
    if root.resolve() != ROOT.resolve() or r.read(root/'CAMPAIGN.json')['campaign'] != CAMPAIGN:
        raise PermissionError('Cleanup is restricted to the approved new campaign root')
    meta = verify_archive(root, tile, draws)
    base = root/'tiles'/tile
    removals = []
    for name in sorted(meta['members']):
        if name.startswith('prepared/'):
            continue
        path = safe_destination(base, name)
        if path.exists():
            if not path.is_file() or r.sha(path) != meta['members'][name]['sha256']:
                raise ValueError('Redundant staging differs; no cleanup performed')
            removals.append(path)
    # Validate the entire removal list first. Never remove an unverified file.
    for path in removals:
        path.unlink()
    directories = {p.parent for p in removals}
    for directory in sorted(directories, key=str):
        if not any(directory.iterdir()):
            directory.rmdir()
    for model in {p.parent.parent for p in removals}:
        if model.exists() and not any(model.iterdir()):
            model.rmdir()
    receipt = dict(passed=True, tile=tile, draws=draws, removed_files=len(removals),
                   archive_sha256=meta['archive_sha256'], recoverable_from=str(archive_path(root, tile, draws)),
                   prepared_preserved=True, job=os.environ.get('SLURM_JOB_ID'))
    r.atomic_json(root/'lifecycle'/f'{tile}_{draws}_{time.time_ns()}.json', receipt)
    return receipt


def campaign_sources():
    names = ('loa_survey_runner_20261016.py', 'pack_loa_survey_tile_20261016.py',
             'loa_survey_production_20261016.py', 'test_loa_survey_production_20261016.py',
             'run_loa_survey_production_20261016.sh', 'loa_survey_preview_20261016.slurm',
             'loa_survey_extension_20261016.slurm')
    return {str(HERE/name): r.sha(HERE/name) for name in names}


def check_campaign(root):
    meta = r.read(root/'CAMPAIGN.json')
    if root.resolve() != ROOT.resolve() or meta['campaign'] != CAMPAIGN:
        raise ValueError('Wrong campaign root')
    if meta['sources'] != campaign_sources() or meta['index_sha256'] != r.sha(root/'INDEX.json'):
        raise ValueError('Campaign source/index drift')
    if r.read(root/'INDEX.json')['sources'] != r.source_binding():
        raise ValueError('Scientific implementation drift')
    return meta


def initialize(root):
    r.require_compute()
    if root != ROOT or root.exists():
        raise FileExistsError('Campaign root must be new and exactly registered')
    index = r.index_catalogue(PREFLIGHT)
    if len(index['cases']) != 2472:
        raise ValueError('Unexpected tile count')
    root.mkdir()
    for folder in ('tiles', 'locks', 'packed', 'logs', 'lifecycle', 'shards'):
        (root/folder).mkdir()
    # Copy modest metadata; reuse the immutable checksum-bound catalogue in place.
    shutil.copyfile(PREFLIGHT/'INDEX.json', root/'INDEX.json')
    (root/'catalogue.h5').symlink_to(PREFLIGHT/'catalogue.h5')
    manifest = dict(campaign=CAMPAIGN, sources=campaign_sources(), index_sha256=r.sha(root/'INDEX.json'),
                    phase='four_draw_preview_only', tiles=2472, shards=103, tiles_per_shard=24,
                    max_concurrent_gpus=8, gpu_hour_ceiling=1236, cpu_node_hour_ceiling=10,
                    preview_gpu_hours_ceiling=206, extension_gpu_hours_ceiling=1030,
                    cpu_limits=dict(prepare_and_validation=2, preview_review=1, final_merge=7),
                    retry_policy='no automatic retry or requeue; further attempts need explicit budget accounting',
                    extension_policy='not submitted until complete preview review is recorded',
                    cleanup_policy='only checksum-matched redundant draw files in this root; archives and prepared inputs retained',
                    created_job=os.environ['SLURM_JOB_ID'])
    r.atomic_json(root/'CAMPAIGN.json', manifest)
    packing.snapshot_sources(root)
    # Separately preserve the production lifecycle/launch code, including untracked files.
    with tarfile.open(root/'PRODUCTION_SOURCES.tar', 'x') as archive:
        for filename in manifest['sources']:
            archive.add(filename, arcname=Path(filename).name, recursive=False)
    r.atomic_json(root/'PRODUCTION_SOURCES.json', dict(sources=manifest['sources'],
                  sha256=r.sha(root/'PRODUCTION_SOURCES.tar')))
    return manifest


def prepare_all(root):
    r.require_compute(); check_campaign(root)
    index = r.read(root/'INDEX.json')
    def worker(number):
        tiles = ','.join(x['id'] for x in index['cases'][number::8])
        command = [PYTHON, str(HERE/'loa_survey_runner_20261016.py'), 'prepare', '--root', str(root),
                   '--tiles', tiles, '--seconds', '6600']
        with (root/'logs'/f'prepare_{number}.log').open('x') as output:
            subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=True)
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(worker, range(8)))
    records = [r.verified_bundle(root/'tiles'/row['id']/'prepared', r.digest(index)) for row in index['cases']]
    if not all(records):
        raise ValueError('Incomplete preparation')
    r.atomic_json(root/'PREPARATION_COMPLETE.json', dict(passed=True, tiles=len(records),
        retained=sum(x['case']['retained_rows'] for x in records),
        property_eligible=sum(x['property_eligible'] for x in records),
        job=os.environ['SLURM_JOB_ID'], index_sha256=r.sha(root/'INDEX.json')))


def validate_recovery(root):
    """Round-trip the preflight's real 32-draw archive without changing it."""
    r.require_compute(); check_campaign(root)
    tile = 'NGC_005_008_003'
    testroot = root/'RECOVERY_VALIDATION'
    testroot.mkdir()
    shutil.copyfile(root/'INDEX.json', testroot/'INDEX.json')
    target = archive_path(testroot, tile, 32)
    target.parent.mkdir()
    shutil.copytree(archive_path(PREFLIGHT, tile, 32), target)
    meta = restore(testroot, tile, 32)
    before = {name: (testroot/'tiles'/tile/name).stat().st_mtime_ns for name in meta['members']}
    restore(testroot, tile, 32)
    for name, value in meta['members'].items():
        path = testroot/'tiles'/tile/name
        if r.sha(path) != value['sha256'] or path.stat().st_mtime_ns != before[name]:
            raise ValueError('Archive restore/no-op replay mismatch')
    r.atomic_json(root/'RECOVERY_VALIDATED.json', dict(passed=True, tile=tile,
        draws=32, members=len(meta['members']), exact_hash_roundtrip=True,
        noop_preserves_mtimes=True, original_preflight_preserved=True,
        archive_sha256=meta['archive_sha256'], job=os.environ['SLURM_JOB_ID']))


def shard(root, number, draws):
    r.require_compute(); check_campaign(root)
    if not 0 <= number < 103 or draws not in (4, 32):
        raise ValueError('Unregistered shard/ensemble')
    if not r.read(root/'PREPARATION_COMPLETE.json')['passed']:
        raise ValueError('Preparation incomplete')
    if not r.read(root/'RECOVERY_VALIDATED.json')['passed']:
        raise ValueError('Recovery validation incomplete')
    if draws == 32:
        review = r.read(root/'PREVIEW_REVIEW.json')
        if not review.get('extension_authorized') or review['campaign_sha256'] != r.sha(root/'CAMPAIGN.json'):
            raise PermissionError('Full preview review must precede extension')
    index = r.read(root/'INDEX.json')
    cases = index['cases'][number*24:(number+1)*24]
    pending = []
    with r.lock(root/'locks'/f'shard_{number}.lock'):
        for row in cases:
            tile = row['id']
            if archive_path(root, tile, draws).exists():
                with r.lock(root/'locks'/f'{tile}.lock'):
                    cleanup(root, tile, draws)
                continue
            if draws == 32:
                with r.lock(root/'locks'/f'{tile}.lock'):
                    restore(root, tile, 4)
            pending.append(tile)
        if pending:
            command = [PYTHON, str(HERE/'loa_survey_runner_20261016.py'), 'infer', '--root', str(root),
                       '--tiles', ','.join(pending), '--first', '0' if draws == 4 else '4', '--stop', str(draws),
                       '--seconds', '6600' if draws == 4 else '35400']
            subprocess.run(command, check=True)
            for tile in pending:
                packing.pack(root, tile, draws)
                with r.lock(root/'locks'/f'{tile}.lock'):
                    cleanup(root, tile, draws)
        check_campaign(root)
        r.atomic_json(root/'shards'/f'{number:03d}_{draws:02d}.json', dict(passed=True,
            tiles=[x['id'] for x in cases], draws=draws, job=os.environ['SLURM_JOB_ID'],
            campaign_sha256=r.sha(root/'CAMPAIGN.json')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['initialize', 'validate-recovery', 'prepare', 'shard'])
    parser.add_argument('--shard', type=int)
    parser.add_argument('--draws', type=int, choices=[4, 32], default=4)
    args = parser.parse_args()
    if args.action == 'initialize': initialize(ROOT)
    elif args.action == 'validate-recovery': validate_recovery(ROOT)
    elif args.action == 'prepare': prepare_all(ROOT)
    else: shard(ROOT, args.shard, args.draws)


if __name__ == '__main__':
    main()
