"""Scientific invariants for release/footprint control diagnostics."""
from pathlib import Path
import sys, unittest
import numpy as np
import healpy as hp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'workflows/p12a_vac'))
from alignment_dr1_controls import quality,count
from alignment_boundary_controls import erode

class Controls(unittest.TestCase):
    def test_quality_half_open_range_and_joint_flags(self):
        a=np.array([(.15,1,1,1),(.55,1,1,1),(.3,0,1,1),(.3,1,0,1),(.3,1,1,0)],dtype=[('z','f8'),('zw','?'),('gal','?'),('dc','?')])
        np.testing.assert_array_equal(quality(a),[1,0,0,0,0])
    def test_nested_pixel_parent_preserves_direction(self):
        ra=np.array([0,22,179,281]);dec=np.array([-40,0,32,78])
        np.testing.assert_array_equal(hp.ang2pix(512,ra,dec,lonlat=True,nest=True)//4,hp.ang2pix(256,ra,dec,lonlat=True,nest=True))
    def test_no_artificial_erosion_at_healpix_vertices(self):
        a=np.ones(hp.nside2npix(256),bool)
        self.assertTrue(erode(a).all())
    def test_isolated_pixel_removed(self):
        a=np.zeros(hp.nside2npix(256),bool);a[12345]=True
        self.assertFalse(erode(a).any())
    def test_accepted_redshift_migration_changes_bins_not_total(self):
        a=np.array([(.20,0),(.34,1)],dtype=[('z','f8'),('cap','i1')]);b=a.copy();b['z']=[.21,.35]
        delta=count(b,np.ones(2,bool))-count(a,np.ones(2,bool))
        self.assertEqual(delta.sum(),0);self.assertGreater(np.abs(delta).sum(),0)
if __name__=='__main__':unittest.main()
