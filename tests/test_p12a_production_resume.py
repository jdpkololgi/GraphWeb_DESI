import importlib.util,sys,tempfile,unittest,json,hashlib
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows/p12a_vac'))
import production_inference as p

class ResumeContract(unittest.TestCase):
    def fixture(self,root):
        (root/'shards').mkdir();case=dict(core_id=23,cap=0,core_start=[0,0,0]);ready=dict(cases=[case]);p.save(root/'INPUTS_READY.json',ready)
        for name in ['a.fits','a.npz']:(root/'shards'/name).write_bytes(b'complete')
        marker=root/'shards/core_000023_COMPLETE.json';record=dict(core_id=23,case=case,input_sha256=p.digest(root/'INPUTS_READY.json'),fits='shards/a.fits',draws='shards/a.npz',fits_sha256=p.digest(root/'shards/a.fits'),draws_sha256=p.digest(root/'shards/a.npz'));p.save(marker,record);return marker,ready
    def test_resume_verifies_completed_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            marker,ready=self.fixture(Path(tmp));self.assertEqual(p.verify_core(marker,ready)['core_id'],23)
            (Path(tmp)/'shards/a.npz').write_bytes(b'truncated')
            with self.assertRaisesRegex(ValueError,'corrupt'):p.verify_core(marker,ready)
    def test_changed_input_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);marker,ready=self.fixture(root);p.save(root/'INPUTS_READY.json',dict(cases=[]))
            with self.assertRaisesRegex(ValueError,'input changed'):p.verify_core(marker,ready)
    def test_changed_ownership_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            marker,ready=self.fixture(Path(tmp));ready['cases'][0]['cap']=1
            with self.assertRaisesRegex(ValueError,'ownership changed'):p.verify_core(marker,ready)
if __name__=='__main__':unittest.main()
