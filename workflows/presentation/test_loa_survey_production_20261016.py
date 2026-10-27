"""Small archive fixtures exercise lifecycle failures; no production data touched."""
import hashlib
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import loa_survey_production_20261016 as p
r = p.r


def fixture(root, draws=4):
    tile = 'TEST_TILE'
    index = dict(cases=[dict(id=tile)], checkpoints={})
    for seed in r.SEEDS:
        index['checkpoints'][f'coarse_{seed}'] = dict(sha256=f'coarse{seed}')
        index['checkpoints'][f'fine_{seed}'] = dict(sha256=f'fine{seed}')
        for a, v in r.FAMILIES[1:]:
            index['checkpoints'][f'{a}_{v}_{seed}'] = dict(sha256=f'{a}{seed}')
    r.atomic_json(root/'INDEX.json', index)
    r.atomic_json(root/'CAMPAIGN.json', dict(campaign=p.CAMPAIGN))
    (root/'lifecycle').mkdir()
    contents = {}
    def bundle(name, receipt):
        payload = ('fake test payload '+name).encode()
        contents[name+'/arrays.npz'] = payload
        receipt.update(tile=tile, passed=True, sha256=hashlib.sha256(payload).hexdigest())
        contents[name+'/COMPLETE.json'] = json.dumps(receipt).encode()
        return receipt
    bundle('prepared', dict(binding=r.digest(index)))
    for seed in r.SEEDS:
        for draw in range(draws):
            coarse = bundle(f'coarse_{seed}/{draw:03d}', dict(binding=f'coarse-{seed}',
                            seed=seed, draw=draw, checkpoint_sha256=f'coarse{seed}'))
            for a, v in r.FAMILIES:
                model = f'{a}_{v}_{seed}'
                bundle(f'{model}/{draw:03d}', dict(seed=seed, draw=draw, model=model,
                    binding=r.digest([coarse['binding'], model, draw, coarse['sha256']]),
                    coarse_sha256=coarse['sha256'], truth_opened=False,
                    checkpoint_sha256=f'fine{seed}' if a == 'legacy' else f'{a}{seed}'))
    destination = p.archive_path(root, tile, draws)
    destination.mkdir(parents=True)
    with tarfile.open(destination/'tile.tar', 'w') as archive:
        for name, data in contents.items():
            info = tarfile.TarInfo(name); info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    meta = dict(passed=True, tile=tile, draws_per_fit=draws, index_sha256=r.sha(root/'INDEX.json'),
                archive_sha256=r.sha(destination/'tile.tar'), members={
                    name:dict(bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
                    for name, data in contents.items()})
    r.atomic_json(destination/'COMPLETE.json', meta)
    return tile, meta


class ProductionTests(unittest.TestCase):
    def test_expected_counts(self):
        self.assertEqual(len(p.expected_members(4)), 66)
        self.assertEqual(len(p.expected_members(32)), 514)
        with self.assertRaises(ValueError): p.expected_members(3)

    def test_realistic_roundtrip_noop_cleanup_and_recovery(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); tile, meta = fixture(root)
            p.restore(root, tile, 4)
            before = {n:(root/'tiles'/tile/n).stat().st_mtime_ns for n in meta['members']}
            p.restore(root, tile, 4)
            self.assertEqual(before, {n:(root/'tiles'/tile/n).stat().st_mtime_ns for n in meta['members']})
            with patch.object(p, 'ROOT', root):
                result = p.cleanup(root, tile, 4)
            self.assertEqual(result['removed_files'], 64)
            self.assertTrue((root/'tiles'/tile/'prepared/arrays.npz').exists())
            p.restore(root, tile, 4)
            for name, row in meta['members'].items():
                self.assertEqual(r.sha(root/'tiles'/tile/name), row['sha256'])

    def test_cleanup_wrong_root_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile, _=fixture(root)
            with self.assertRaises(PermissionError): p.cleanup(root,tile,4)

    def test_cleanup_prevalidates_all_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile, _=fixture(root); p.restore(root,tile,4)
            bad=root/'tiles'/tile/'unet_x_29/003/arrays.npz'
            bad.write_bytes(b'corrupt fixture')
            with patch.object(p,'ROOT',root), self.assertRaises(ValueError): p.cleanup(root,tile,4)
            self.assertTrue((root/'tiles'/tile/'coarse_17/000/arrays.npz').exists())

    def test_partial_cleanup_recoverable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile, meta=fixture(root); p.restore(root,tile,4)
            missing=root/'tiles'/tile/'coarse_17/000/arrays.npz'; missing.unlink()
            p.restore(root,tile,4)
            self.assertEqual(r.sha(missing),meta['members']['coarse_17/000/arrays.npz']['sha256'])

    def test_existing_corruption_preserved_and_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile,_=fixture(root); p.restore(root,tile,4)
            bad=root/'tiles'/tile/'prepared/arrays.npz'; bad.write_bytes(b'corrupt fixture')
            with self.assertRaises(ValueError): p.restore(root,tile,4)
            self.assertEqual(bad.read_bytes(),b'corrupt fixture')

    def test_archive_corruption_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile,_=fixture(root)
            with (p.archive_path(root,tile,4)/'tile.tar').open('ab') as stream: stream.write(b'corruption')
            with self.assertRaises(ValueError): p.restore(root,tile,4)
            self.assertFalse((root/'tiles').exists())

    def test_unsafe_paths_and_symlinks(self):
        with tempfile.TemporaryDirectory() as temporary:
            base=Path(temporary)/'staging'; base.mkdir()
            for name in ('../outside','/tmp/outside','a/../../outside'):
                with self.assertRaises(ValueError): p.safe_destination(base,name)
            (base/'link').symlink_to(Path(temporary))
            with self.assertRaises(ValueError): p.safe_destination(base,'link/file')

    def test_archive_membership_mismatch_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile,meta=fixture(root)
            meta['members']['../escape'] = dict(bytes=0,sha256='bad')
            (p.archive_path(root,tile,4)/'COMPLETE.json').write_text(json.dumps(meta))
            with self.assertRaises(ValueError): p.restore(root,tile,4)

    def test_index_mismatch_refused(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary); tile,_=fixture(root)
            (root/'INDEX.json').write_text('{}')
            with self.assertRaises(ValueError): p.restore(root,tile,4)


if __name__ == '__main__':
    unittest.main(verbosity=2)
