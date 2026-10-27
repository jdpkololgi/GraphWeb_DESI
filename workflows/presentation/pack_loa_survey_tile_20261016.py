"""Export a complete tile to one verified archive; never delete draw bundles.

For future production, use job-local staging plus an explicitly approved
cleanup policy to avoid retaining millions of small filesystem entries. This
export alone does not authorize deletion or change the inference resume path.
Archives contain private galaxy data and must not be embedded in public Pages.
"""
import argparse
import hashlib
import os
from pathlib import Path
import tarfile
import tempfile
import time
import loa_survey_runner_20261016 as runner


def snapshot_sources(root):
    """Keep the exact inspected code beside the receipts, not hashes alone."""
    index=runner.read(root/'INDEX.json'); sources=dict(index['sources'])
    for path, expected in sources.items():
        if runner.sha(path)!=expected:
            raise ValueError('Source changed before snapshot: '+path)
    for name in ('e2e_coupled_data_v1.json','e2e_coupled_coordinates_v2.json','e2e_coupled_product_layout_v1.json'):
        path=runner.ILLUSTRIS/'configs'/name
        sources[str(path)]=runner.sha(path)
    for name in ('loa_survey_preflight_checks_20261016.py','loa_survey_empty_tile_check_20261016.py',
                 'test_loa_survey_runner_20261016.py','run_loa_survey_preflight_20261016.sh',
                 'pack_loa_survey_tile_20261016.py'):
        path=Path(__file__).with_name(name)
        sources[str(path)]=runner.sha(path)
    output=root/'SOURCE_SNAPSHOT'
    if output.exists():
        previous=runner.read(output/'COMPLETE.json')
        if not previous['passed'] or previous['sources']!=sources or runner.sha(output/'sources.tar')!=previous['sha256']:
            raise ValueError('Existing source snapshot differs')
        return
    temporary=Path(tempfile.mkdtemp(prefix='.SOURCE_SNAPSHOT.',dir=root))
    with tarfile.open(temporary/'sources.tar','w') as archive:
        for filename in sorted(sources):
            path=Path(filename); info=tarfile.TarInfo(filename.lstrip('/'))
            info.size=path.stat().st_size; info.mode=0o640
            with path.open('rb') as stream:
                archive.addfile(info,stream)
    with tarfile.open(temporary/'sources.tar','r') as archive:
        for member in archive:
            with archive.extractfile(member) as stream:
                if hashlib.sha256(stream.read()).hexdigest()!=sources['/'+member.name]:
                    raise ValueError('Source archive checksum mismatch')
    runner.atomic_json(temporary/'COMPLETE.json',dict(passed=True,sources=sources,
        sha256=runner.sha(temporary/'sources.tar'),job=os.environ['SLURM_JOB_ID'],
        note='Source hashes checked against INDEX; additional launcher/test/config files recorded at closeout'))
    os.rename(temporary,output)


def pack(root, tile, draws):
    runner.require_compute()
    start=time.monotonic(); index=runner.read(root/'INDEX.json')
    if tile not in {r['id'] for r in index['cases']}:
        raise ValueError('Unknown tile')
    source=runner.digest(index); base=root/'tiles'/tile
    output=root/'packed'/f'{tile}_{draws:02d}draws'
    output.parent.mkdir(parents=True,exist_ok=True)
    with runner.lock(root/'locks'/(tile+'.lock')):
        if output.exists():
            raise FileExistsError('Archive exists; no overwrite: '+str(output))
        runner.verified_bundle(base/'prepared',source)
        bundles=[base/'prepared']
        for seed in runner.SEEDS:
            for draw in range(draws):
                coarse=base/f'coarse_{seed}'/f'{draw:03d}'
                record=runner.read(coarse/'COMPLETE.json')
                runner.verified_bundle(coarse,record['binding'])
                bundles.append(coarse)
                for architecture,prediction in runner.FAMILIES:
                    model=f'{architecture}_{prediction}_{seed}'
                    path=base/model/f'{draw:03d}'
                    receipt=runner.verified_bundle(path,runner.digest([record['binding'],model,draw,record['sha256']]))
                    if receipt is None:
                        raise ValueError('Missing draw: '+str(path))
                    bundles.append(path)
        files=[p/name for p in bundles for name in ('arrays.npz','COMPLETE.json')]
        temporary=Path(tempfile.mkdtemp(prefix='.'+output.name+'.',dir=output.parent))
        members={}
        with tarfile.open(temporary/'tile.tar',mode='w',format=tarfile.PAX_FORMAT) as archive:
            for path in files:
                name=str(path.relative_to(base))
                info=tarfile.TarInfo(name); info.size=path.stat().st_size; info.mode=0o640
                with path.open('rb') as stream:
                    archive.addfile(info,stream)
                members[name]={'bytes':path.stat().st_size,'sha256':runner.sha(path)}
        # Verify every payload after serialization, not only the container hash.
        with tarfile.open(temporary/'tile.tar',mode='r') as archive:
            if set(archive.getnames()) != set(members) or len(archive.getmembers()) != len(members):
                raise ValueError('Archive member list mismatch')
            for member in archive:
                if not member.isfile():
                    raise ValueError('Non-file member in archive')
                h=hashlib.sha256()
                with archive.extractfile(member) as stream:
                    for block in iter(lambda:stream.read(8<<20),b''):
                        h.update(block)
                if h.hexdigest()!=members[member.name]['sha256']:
                    raise ValueError('Archive payload checksum mismatch')
        result=dict(passed=True,tile=tile,draws_per_fit=draws,index_sha256=runner.sha(root/'INDEX.json'),
                    archive_sha256=runner.sha(temporary/'tile.tar'),bytes=(temporary/'tile.tar').stat().st_size,
                    members=members,seconds=time.monotonic()-start,job=os.environ['SLURM_JOB_ID'],
                    source_sha256=runner.sha(__file__),input_bundles_preserved=True,
                    contains_private_galaxy_data=True)
        runner.atomic_json(temporary/'COMPLETE.json',result)
        # Flush both payload and directory before publishing an archive that may
        # subsequently become the authoritative copy for production recovery.
        with (temporary/'tile.tar').open('rb') as stream:
            os.fsync(stream.fileno())
        descriptor=os.open(temporary,os.O_RDONLY|os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        os.rename(temporary,output)
        descriptor=os.open(output.parent,os.O_RDONLY|os.O_DIRECTORY)
        try:
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        print({k:v for k,v in result.items() if k!='members'},flush=True)
        return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=runner.DEFAULT_ROOT)
    parser.add_argument('--tile',required=True)
    parser.add_argument('--draws',type=int,choices=[4,32],required=True)
    args=parser.parse_args()
    pack(args.root,args.tile,args.draws)
    snapshot_sources(args.root)
