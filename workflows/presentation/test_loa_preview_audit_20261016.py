"""Small numerical validation tests; no survey input or GPU required."""
import unittest
import numpy as np
from audit_loa_survey_preview_20261016 import check_arrays


def arrays(rows=2):
    return dict(rho=np.ones((32,16,16)),tensor=np.zeros((32,16,16,6)),
                galaxy_eigen=np.zeros((rows,3)),galaxy_rho=np.ones(rows))


class AuditTests(unittest.TestCase):
    def test_valid(self): check_arrays(arrays(),2)
    def test_empty_galaxies(self): check_arrays(arrays(0),0)
    def test_nonfinite(self):
        value=arrays(); value['rho'][0,0,0]=np.nan
        with self.assertRaises(ValueError): check_arrays(value,2)
    def test_density(self):
        value=arrays(); value['rho'][0,0,0]=0
        with self.assertRaises(ValueError): check_arrays(value,2)
    def test_order(self):
        value=arrays(); value['galaxy_eigen'][0]=[1,0,-1]
        with self.assertRaises(ValueError): check_arrays(value,2)
    def test_trace(self):
        value=arrays(); value['tensor'][0,0,0,0]=1
        with self.assertRaises(AssertionError): check_arrays(value,2)
    def test_shape(self):
        with self.assertRaises(ValueError): check_arrays(arrays(1),2)
    def test_mass(self): check_arrays(arrays(),2,np.ones((8,4,4)))
    def test_mass_mismatch(self):
        with self.assertRaises(AssertionError): check_arrays(arrays(),2,np.full((8,4,4),2.))


if __name__=='__main__': unittest.main(verbosity=2)
