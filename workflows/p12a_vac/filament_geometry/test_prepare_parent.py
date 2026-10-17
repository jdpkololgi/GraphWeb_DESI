"""Small geometry contracts; no catalogue scan or scheduler calls."""
import unittest
import numpy as np
from prepare_parent import coordinates,rotation,CENTER,RADIUS
class Geometry(unittest.TestCase):
 def test_axes(self):
  np.testing.assert_allclose(coordinates(np.array([180.,270.,180.]),np.array([0.,0.,90.]),np.ones(3)),[[0,1,0],[1,0,0],[0,0,1]],atol=1e-15)
 def test_rotation_and_padding(self):
  corner=np.array([180.,225.,45.]);r=1.5*7/.6766
  for angle in [0,30,60,90]:
   q=rotation(angle);np.testing.assert_allclose(q@q.T,np.eye(3),atol=1e-15)
   self.assertLess(np.linalg.norm(corner@q)+r,RADIUS)
 def test_unrotated_selection(self):
  p=np.array([[0,1025,0],[0,1025,11.]])-CENTER
  q=p@rotation(0).T
  np.testing.assert_array_equal(abs(q[:,2])<10,[True,False])
if __name__=='__main__':unittest.main()
