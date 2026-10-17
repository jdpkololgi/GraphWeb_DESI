import unittest
import numpy as np
from compare_offslab import shape,contrast2
class ShapeTests(unittest.TestCase):
 def test_shapes(self):
  line=np.array([[-1,0,0],[0,0,0],[1,0,0]])
  plane=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]])
  ball=np.r_[plane,[[0,0,1],[0,0,-1]]]
  for x,expected in [(line,[1,0,0]),(plane,[0,1,0]),(ball,[0,0,1])]:np.testing.assert_allclose(shape(x)[0],expected)
 def test_rotation(self):
  rng=np.random.default_rng(5);x=rng.normal(size=(30,3))*[3,2,.3];q,_=np.linalg.qr(rng.normal(size=(3,3)))
  np.testing.assert_allclose(shape(x)[0],shape(x@q+19)[0],atol=1e-12)
 def test_line_slice_of_wall(self):
  g=np.arange(-5,6);x,y=np.meshgrid(g,g);wall=np.c_[x.ravel(),np.zeros(x.size),y.ravel()]
  self.assertGreater(contrast2(wall[abs(wall[:,2])<.5,:2]),.99)
  self.assertGreater(shape(wall)[0][1],.99)
if __name__=='__main__':unittest.main()
