import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.PolyFEMApplication import *
kernel = Kernel()   #defining kernel

util = PolyTreeUtility()
util.TestPolygonDecomposition(100, 1.0e-4, 3)

