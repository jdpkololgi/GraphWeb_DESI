import unittest,importlib.util
from pathlib import Path
import numpy as np
s=importlib.util.spec_from_file_location('observer',Path(__file__).resolve().parents[1]/'workflows/catalog/p12a_observation_patch.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
class PatchTests(unittest.TestCase):
    def test_cic_conservation(self):
        x=np.array([[1.,1.,1.],[1.25,1.5,1.75]])
        a=m.count_patch(x,np.array([0,1]),np.zeros(3),1.,np.zeros(3,dtype=int),np.array([4,4,4]))
        self.assertAlmostEqual(float(a.sum()),2.)
    def test_partition_keeps_neighbor_context(self):
        x=np.array([[1.8,1.,1.],[2.2,1.,1.]])
        whole=m.count_patch(x,np.array([0,0]),np.zeros(3),1.,np.zeros(3,dtype=int),np.array([4,4,4]))
        left=m.count_patch(x,np.array([0,0]),np.zeros(3),1.,np.zeros(3,dtype=int),np.array([2,4,4]))
        right=m.count_patch(x,np.array([0,0]),np.zeros(3),1.,np.array([2,0,0]),np.array([4,4,4]))
        np.testing.assert_array_equal(np.concatenate([left,right]),whole)
    def test_success_policy(self):
        d=np.array([(0.2,0,25.,'GALAXY'),(0.2,0,24.,'GALAXY'),(0.2,0,50.,'STAR'),(float('nan'),0,50.,'GALAXY')],dtype=[('Z_not4clus','f8'),('ZWARN','i8'),('DELTACHI2','f8'),('SPECTYPE','U8')])
        np.testing.assert_array_equal(m.successful_rows(d),[True,False,False,False])
        np.testing.assert_array_equal(m.successful_rows(d,mock=True),[True,True,True,False])
    def test_invalid_coordinates(self):
        with self.assertRaises(ValueError):m.observer_xyz([0.],[0.],[float('nan')])
        with self.assertRaises(ValueError):m.observer_xyz([0.],[0.],[.9])
if __name__=='__main__':unittest.main()
