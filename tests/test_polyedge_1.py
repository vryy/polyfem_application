import sys
import os
kratos_root_path=os.environ['KRATOS_ROOT_PATH']
##################################################################
##################################################################
#importing Kratos modules
from KratosMultiphysics import *
from KratosMultiphysics.PolyFEMApplication import *
kernel = Kernel()   #defining kernel

node1 = PolyVertex2D(1, 1.0, 0.0)
node2 = PolyVertex2D(2, 1.0, 1.0)

edge = PolyHalfEdge2D(node1, node2)
print(str(edge))
print(edge.HashCode())

