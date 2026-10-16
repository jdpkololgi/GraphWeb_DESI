"""Small arithmetic tests; no catalogue reads or scheduler actions."""
import unittest
import numpy as np
from run import limiting_curve, selection, controlled, unique_only


class SelectionTests(unittest.TestCase):
    def test_mass_rescaling_and_conservative_population(self):
        n=2000
        z=np.full(n,.205);m=np.r_[np.full(n//2,10.),np.full(n//2,11.)]
        r=np.full(n,19.);lim=np.full(n,19.5)
        sf=np.r_[np.full(n//2,10.),np.zeros(n//2)]
        curve=limiting_curve(z,m,r,lim,sf,np.array([.2,.21]))
        self.assertAlmostEqual(curve[0]['logmass_limit'],10.8)
        self.assertEqual(curve[0]['populations']['quiescent']['rows'],1000)

    def test_sparse_population_fails_closed(self):
        c=limiting_curve(np.full(101,.205),np.full(101,10.),np.full(101,19.),np.full(101,19.5),np.ones(101),np.array([.2,.21]))
        self.assertFalse(c[0]['usable'])
        self.assertIsNone(c[0]['logmass_limit'])

    def test_fixed_threshold_and_boundaries(self):
        curve=[{'zlo':round(.2+i*.01,2),'zhi':round(.21+i*.01,2),'usable':True,'logmass_limit':10.5+i*.01} for i in range(10)]
        z=np.array([.199,.2,.299,.3]);m=np.full(4,10.60)
        _,vary,fixed,cut=selection(z,m,curve)
        np.testing.assert_array_equal(fixed,[False,True,True,False])
        np.testing.assert_array_equal(vary,fixed)
        self.assertAlmostEqual(cut,10.59)

    def test_duplicate_ids_not_arbitrarily_chosen(self):
        a=np.array([(1,),(1,),(2,)],dtype=[('TARGETID','i8')]);b,n=unique_only(a)
        self.assertEqual(n,2);self.assertEqual(b['TARGETID'].tolist(),[2])

    def test_common_mix_recovers_known_contrast(self):
        group=np.tile(np.arange(4),400);cell=np.repeat([0,1],800)
        y=group+10*cell;w=np.eye(4)[group]
        r=controlled(y,{'a':w,'b':w},cell,np.arange(1600)%17,nboot=32)
        self.assertEqual(r['status'],'measured')
        self.assertAlmostEqual(r['methods']['a']['contrast'],3.)
        self.assertEqual(r['methods']['a'],r['methods']['b'])

    def test_no_common_support(self):
        n=40;w=np.zeros((n,4));w[:,0]=1
        r=controlled(np.ones(n),{'a':w},np.zeros(n),np.arange(n)%4)
        self.assertIn('unavailable',r['status'])


if __name__=='__main__':unittest.main()
