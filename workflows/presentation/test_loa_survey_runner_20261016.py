"""Focused integrity/selection tests; run on the approved CPU allocation."""
import importlib.util
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

spec = importlib.util.spec_from_file_location('survey',Path(__file__).with_name('loa_survey_runner_20261016.py'))
survey = importlib.util.module_from_spec(spec)
spec.loader.exec_module(survey)


class RunnerTests(unittest.TestCase):
    def test_seed_addresses(self):
        values = [survey.draw_seed(t,s,d,k) for t in ('NGC_0','SGC_0')
                  for s in survey.SEEDS for d in range(32) for k in ('coarse','fine','legacy-fine')]
        self.assertEqual(len(values),len(set(values)))
        self.assertEqual(survey.draw_seed('tile',17,3,'fine'),survey.draw_seed('tile',17,3,'fine'))

    def test_shell_edges(self):
        self.assertEqual([survey.shell_index(z) for z in (.15,.249,.25,.35,.45,.55)], [0,0,1,2,3,3])
        with self.assertRaises(ValueError):
            survey.shell_index(.551)

    def test_geometry_selection(self):
        cases = [dict(id=f'{cap}_{shell}_{i}',cap=cap,shell=shell,retained_rows=1,boundary_median=i)
                 for cap in ('NGC','SGC') for shell in range(4) for i in range(9)]
        selected = survey.select_preflight(cases)
        self.assertEqual(len(selected),16)
        self.assertTrue(all(r['boundary_median']==(2 if r['stratum']=='edge' else 6) for r in selected))
        self.assertEqual(survey.select_preflight(list(reversed(cases))),selected)

    def test_shards_exactly_once(self):
        cases = [dict(id=str(i)) for i in range(2472)]
        found=[]
        for shard in range(103):
            args=SimpleNamespace(scope='all',tiles=None,cap=None,shard=shard)
            found.extend(survey.selected_cases(dict(cases=cases),args))
        self.assertEqual(found,cases)
        args.shard=103
        with self.assertRaises(ValueError):
            survey.selected_cases(dict(cases=cases),args)

    def test_atomic_bundle_and_resume(self):
        import numpy as np
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'draw'
            self.assertIsNone(survey.verified_bundle(path,'binding'))
            receipt=survey.commit_bundle(path,{'a':np.arange(5)},dict(draw=0),'binding')
            self.assertEqual(receipt,survey.verified_bundle(path,'binding'))
            with self.assertRaises(ValueError):
                survey.verified_bundle(path,'changed')
            with self.assertRaises(FileExistsError):
                survey.commit_bundle(path,{'a':np.arange(5)},dict(draw=0),'binding')
            with (path/'arrays.npz').open('ab') as f:
                f.write(b'corruption')
            with self.assertRaises(ValueError):
                survey.verified_bundle(path,'binding')

    def test_interrupted_temporary_not_complete(self):
        with tempfile.TemporaryDirectory() as temporary:
            parent=Path(temporary)
            (parent/'.draw.incomplete').mkdir()
            self.assertIsNone(survey.verified_bundle(parent/'draw','binding'))
            (parent/'draw').mkdir()
            with self.assertRaises(FileNotFoundError):
                survey.verified_bundle(parent/'draw','binding')

    def test_json_no_clobber(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'receipt.json'
            survey.atomic_json(path,{'value':1})
            with self.assertRaises(FileExistsError):
                survey.atomic_json(path,{'value':2})
            self.assertEqual(json.loads(path.read_text()),{'value':1})

    def test_tile_lock(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'tile.lock'
            with survey.lock(path):
                with self.assertRaises(BlockingIOError):
                    with survey.lock(path):
                        pass
            with survey.lock(path):
                pass

    def test_shared_index_readers(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/'index.lock'
            with survey.lock(path,shared=True):
                with survey.lock(path,shared=True):
                    pass
                with self.assertRaises(BlockingIOError):
                    with survey.lock(path):
                        pass


if __name__=='__main__':
    unittest.main(verbosity=2)
