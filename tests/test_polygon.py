import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.PolyFEMApplication import *
kernel = Kernel()   #defining kernel

util = PolyFEMUtility()
util.TestPolygonShapeFunction(4, 0.3, 0.1)
util.TestPolygonShapeFunction(5, 0.0, 0.0)
util.TestPolygonShapeFunction(5, 0.1, 0.2)
util.TestPolygonShapeFunction(6, 0.5, 0.4)

util.TestPolygonQuadrature(5, 1)
util.TestPolygonQuadrature(5, 2)

