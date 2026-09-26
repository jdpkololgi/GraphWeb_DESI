import unittest,tempfile,importlib.util
from pathlib import Path
import numpy as np
import fitsio
s=importlib.util.spec_from_file_location('trial',Path(__file__).resolve().parents[1]/'workflows/p12a_vac/loa_trial.py');t=importlib.util.module_from_spec(s);s.loader.exec_module(t)
class TrialTests(unittest.TestCase):
 def test_contiguous_scan_matches_canonical(self):
  d=np.array([(130.,20.,'N',True,0),(130.,20.,'N',False,0),(10.,-20.,'S',True,4),(131.,21.,'S',True,0)],dtype=[('RA','f8'),('DEC','f8'),('PHOTSYS','S1'),('GOODHARDLOC','?'),('MASKBITS','i8')])
  with tempfile.TemporaryDirectory() as root:
   p=Path(root)/'random.fits';fitsio.write(p,d)
   a=np.zeros((4,t.hp.nside2npix(t.NSIDE)),dtype='i8');b=np.zeros_like(a)
   x=t.add_random_file(a,{'path':str(p)});y=t.contiguous_random_file(b,{'path':str(p)})
   np.testing.assert_array_equal(a,b)
   for key in ['rows','accepted_rows','goodhardloc_rejected','maskbits_nonzero']:self.assertEqual(x[key],y[key])
if __name__=='__main__':unittest.main()
