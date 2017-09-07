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
util.TestClusterPoints1()
util.TestClusterPoints2()
util.TestClusterLengths1()

