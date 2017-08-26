import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.PolyFEMApplication import *
kernel = Kernel()   #defining kernel

vertex1 = PolyVertex2D(1, 0.0, 1.0)
print(str(vertex1))
