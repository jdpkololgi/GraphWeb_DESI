import importlib.util,sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows/p12a_vac'))
from flag_joint_magnitude import intersections
class IntersectionTests(unittest.TestCase):
 def test_overlapping_cuts_are_counted_once(self):
  p=np.zeros((2,256),dtype='i8');p[:,255]=10;p[:,254]=3;p[:,252]=5
  c=intersections(p)
  np.testing.assert_array_equal(c[:,0],[18,18]);np.testing.assert_array_equal(c[:,1],[10,10]);np.testing.assert_array_equal(c[:,2],[13,13]);np.testing.assert_array_equal(c[:,3],[10,10]);np.testing.assert_array_equal(c[:,255],[10,10])
if __name__=='__main__':unittest.main()
