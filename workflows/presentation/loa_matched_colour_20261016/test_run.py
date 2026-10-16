"""Small contracts: remove known confounding; retain tails; discard unsupported cells."""
import unittest
import numpy as np
from run import common_support, estimate

class MatchingTests(unittest.TestCase):
    def test_known_mass_mixture_is_removed_and_tails_retained(self):
        mass=[];z=[];colour=[];classes=[]
        for cell in range(2):
            for env in range(4):
                n=40*(1+(env if cell else 3-env))
                mass.extend([10.65+.2*cell]*n);z.extend([.155]*n)
                colour.extend([-1. if cell==0 else 3.]*n);classes.extend([env]*n)
        m=np.array(mass);z=np.array(z);y=np.array(colour);p=np.eye(4)[classes]
        keep,c,ids,target=common_support(m,z,p,.1,.01)
        self.assertTrue(keep.all())
        out=estimate(y,m,z,p,c,np.zeros(len(y),int),target,np.ones((256,1)))
        self.assertGreater(np.ptp(np.array(out['raw_common'])[:,0]),.1)
        self.assertTrue(np.allclose(out['matched'],out['matched'][0]))
        self.assertAlmostEqual(out['cluster_minus_void'],0)
        self.assertAlmostEqual(sum(out['matched'][0]),1)
        self.assertAlmostEqual(sum(out['matched'][0][1:-1]),0)
        for b in out['balance']:self.assertLess(b['max_cell_error'],1e-12)

    def test_common_support_requires_every_environment(self):
        p=np.eye(4)[np.tile(np.arange(4),30)]
        p=np.vstack((p,np.tile([1,0,0,0],(100,1))))
        m=np.r_[np.full(120,10.65),np.full(100,10.85)]
        keep,c,ids,target=common_support(m,np.full(220,.155),p,.1,.01)
        self.assertEqual(keep.sum(),120)
        self.assertEqual(len(ids),1)
        self.assertEqual(target.tolist(),[1.])

if __name__=='__main__':unittest.main()
