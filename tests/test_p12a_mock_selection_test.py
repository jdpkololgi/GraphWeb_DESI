import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec=importlib.util.spec_from_file_location('selection',Path(__file__).resolve().parents[1]/'workflows/p12a_vac/mock_selection_test.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

class SelectionTests(unittest.TestCase):
    def test_quality_boundaries_are_not_interchanged(self):
        d=np.array([(0,25,'GALAXY'),(0,40,'GALAXY'),(0,41,'STAR'),(4,99,'GALAXY')],dtype=[('ZWARN','i8'),('DELTACHI2','f8'),('SPECTYPE','U8')])
        q=m.quality(d)
        np.testing.assert_array_equal(q['quality25'],[True,True,False,False])
        np.testing.assert_array_equal(q['quality40'],[False,False,True,False])
    def test_brightening_adds_rows_without_duplication(self):
        z=np.array([.2,.5,.5]);r=np.array([19.49,19.6,20.])
        baseline=m.trial_mask(z,r,0,0);trial=m.trial_mask(z,r,0,.4)
        self.assertTrue(np.all(trial[baseline]));np.testing.assert_array_equal(trial,[True,True,False])
    def test_cap_counts_conserve_selected_rows(self):
        z=np.array([.151,.251,.451]);c=np.array([0,1,0]);k=np.array([True,False,True])
        self.assertEqual(m.counts(z,c,k).sum(),2)
        np.testing.assert_array_equal(m.shell_counts(m.counts(z,c,k)),[[1,0,0,1],[0,0,0,0]])
if __name__=='__main__':unittest.main()
