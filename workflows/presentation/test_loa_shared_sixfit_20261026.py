"""Small synthetic geometry/integration tests; no checkpoints or truth."""
import unittest
import numpy as np
import torch
import loa_shared_sixfit_20261026 as s

class SharedTests(unittest.TestCase):
    def test_expanded_geometry(self):
        mids=[[96,96,96],[128,128,128]]
        g=s.geometry(mids,True);c=s.Geometry.make(mids,'coarse')
        self.assertEqual(g.window,(96,64,64))
        for i in range(2):
            fine=np.asarray(g.origin)+g.starts[i]
            coarse=np.asarray(c.origin)+c.starts[i]+[12,16,16]
            np.testing.assert_array_equal(fine,coarse*4)
        self.assertTrue(all(n%4==0 for n in g.shape))

    def test_geometry_guards(self):
        for mid in ([],[[1,8,8]],[[8,8,8],[8,8,8]]):
            with self.assertRaises(ValueError):s.geometry(mid,True)

    def test_consensus_and_uncovered(self):
        class Constant:
            def schedule(self,t):return t
            def __call__(self,x,t,value):return torch.ones_like(x)*value
        g=s.Geometry(((0,0,0),(2,0,0)),(8,4,4),(0,0,0),(4,4,4))
        initial=torch.zeros((1,1,*g.shape))
        result,w=s.integrate(Constant(),'coarse',[1.,3.],initial,g,4)
        self.assertTrue(torch.equal(result[0,0,:,0,0],torch.tensor([1.,1.,2.,2.,3.,3.,0.,0.])))
        self.assertTrue(torch.equal(w[0,0,:,0,0],torch.tensor([1.,1.,2.,2.,1.,1.,0.,0.])))

    def test_fine_block_projection(self):
        class Slope:
            def schedule(self,t):return t
            def __call__(self,x,t,c):return torch.arange(8,dtype=x.dtype).reshape(1,1,8,1,1).expand_as(x)
        g=s.Geometry(((0,0,0),),(8,4,4),(0,0,0),(8,4,4))
        result,_=s.integrate(Slope(),'fine',[None],torch.zeros((1,1,8,4,4)),g,4)
        self.assertLess(float(s.block_mean(result).abs().max()),1e-7)

    def test_expanded_consensus(self):
        class Model:
            def eval(self):return self
            def prepare(self,c):return c
            def velocity(self,x,t,c):return s.project(torch.arange(8,dtype=x.dtype).reshape(1,1,8,1,1).expand_as(x))*c
        g=s.Geometry(((0,0,0),(4,0,0)),(12,4,4),(0,0,0),(8,4,4))
        result,w=s.integrate_expanded(Model(),[1.,3.],torch.zeros((1,1,12,4,4)),g,4)
        self.assertLess(float(s.block_mean(result).abs().max()),1e-7)
        pattern=torch.tensor([-1.5,-.5,.5,1.5])
        torch.testing.assert_close(result[0,0,:,0,0],torch.cat([pattern,pattern*2,pattern*3]))
        bad=s.Geometry(((0,0,0),(2,0,0)),(12,4,4),(0,0,0),(8,4,4))
        with self.assertRaises(ValueError):s.integrate_expanded(Model(),[1.,3.],torch.zeros((1,1,12,4,4)),bad,4)

if __name__=='__main__':unittest.main()
