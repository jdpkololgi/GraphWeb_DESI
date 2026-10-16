"""Small contract checks; no catalogue reads or heavy compute."""
import unittest
import numpy as np
from run import category, index_match, sep_arcsec, raw, controlled

class Contract(unittest.TestCase):
    def test_colour_boundaries_and_invalid(self):
        c=np.array([.749,.75,.751,1.499,1.5,np.nan])
        np.testing.assert_array_equal(category(c),[0,2,1,1,2,2])
        np.testing.assert_array_equal(category(c,.70),[1,1,1,1,2,2])

    def test_match_requires_exact_id(self):
        ids=np.array([10,20,30]);at,ok=index_match(ids,np.array([5,10,25,30,35]))
        np.testing.assert_array_equal(ok,[False,True,False,True,False])
        np.testing.assert_array_equal(ids[at[ok]],[10,30])

    def test_ra_wrap(self):
        self.assertAlmostEqual(float(sep_arcsec(np.array([359.9999]),np.array([0.]),np.array([.0001]),np.array([0.]))[0]),.72,places=6)

    def test_constant_property_has_zero_environment_contrast(self):
        # Every class and sky block is represented; a constant outcome must not
        # acquire a spurious environmental trend through bootstrap or matching.
        p=np.tile(np.eye(4),(80,1));y=np.full(320,2.5)
        region=np.repeat(np.arange(8),40);cells=np.repeat(np.arange(2),160)
        a=raw(y,p,region)
        self.assertEqual(a['bootstrap_valid'],256)
        np.testing.assert_allclose(a['mean'],2.5)
        np.testing.assert_allclose(a['contrast_ci95'],0)
        b=controlled(y,{'soft':p},cells,region)
        self.assertEqual(b['status'],'measured')
        np.testing.assert_allclose(b['methods']['soft']['mean'],2.5)
        np.testing.assert_allclose(b['methods']['soft']['contrast_ci95'],0)

if __name__=='__main__':unittest.main()
