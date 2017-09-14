from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import sys
import math
import re

from KratosMultiphysics import *
from KratosMultiphysics.PolyFEMApplication import *
CheckForPreviousImport()

def IsInt(n):
    try:
        int(n)
        return True
    except ValueError:
        return False

def IsFloat(d):
    try:
        float(d)
        return True
    except ValueError:
        return False

class PolyTreeRefinementUtility:

    def __init__(self, name, model_part):
        self.name = name
        self.tree = PolyTree2D()
        self.tree.Synchronize(model_part, True)
        print("PolyTree " + str(self.name) + " is read successfully from model_part")
        print("################################################")
        print("###### " + self.name + " REFINE APP")
        print("################################################")
        self.help = """Commands:
    h - print this help
    a - probe the parent face (if has) of a specific face, e.g. a,1,2
    r - follow by number to indicate which face will be refined, e.g. r,1,2,3
    c - follow by number to indicate which face will be refined, e.g. c,1,2,3
    f - begin and finalize the refine/coarsen operation
    v - validate the tree
    t - print the tree
    p - export the python file
    m - export the matlab file to visualize the tree
    b - set the MERGE_PARAMETER
    q - quit the application
    Note: all command must be separated by comma"""
        self.operations = []
        self.cmd_history = "h"

    def Loop(self):
        while True:
            sys.stdout.write('/> ')
            cmd = raw_input()
            cmd_str = cmd.split(',')
            # cmd_str = re.split('; |, |\*|\n', cmd)
            # print(cmd_str)
            flag = self.Parse(cmd_str)
            if flag != 0:
                print('Quit!!!')
                break;

    def Parse(self, cmd_str):
        cnt = 0
        while cnt < len(cmd_str):
            dat = cmd_str[cnt]
            if dat == 'q':
                return 1
            if dat == 'h':
                print(self.help)
            elif dat == 'r' or dat == 'c':
                cnt = cnt + 1
                num = []
                while cnt < len(cmd_str):
                    if IsInt(cmd_str[cnt]):
                        num.append(int(cmd_str[cnt]))
                        cnt = cnt + 1
                    else:
                        cnt = cnt - 1
                        break
                if dat == 'r':
                    print("Faces to be refined: " + str(num))
                    self.cmd_history = self.cmd_history + ',r'
                    refined_face = []
                    for f in num:
                        err = self.tree.MarkFaceRefine(f)
                        if err == 0:
                            refined_face.append(f)
                            self.cmd_history = self.cmd_history + ',' + str(f)
                    self.operations.append(['r', refined_face])
                elif dat == 'c':
                    print("Faces to be coarsen: " + str(num))
                    self.cmd_history = self.cmd_history + ',c'
                    coarsen_face = []
                    for f in num:
                        err = self.tree.MarkFaceCoarsen(f)
                        if err == 0:
                            coarsen_face.append(f)
                            self.cmd_history = self.cmd_history + ',' + str(f)
                    self.operations.append(['c', coarsen_face])
            elif dat == 'f':
                self.tree.BeginRefineCoarsen()
                self.tree.EndRefineCoarsen()
                self.cmd_history = self.cmd_history + ",f"
            elif dat == 'b':
                cnt = cnt + 1
                if IsFloat(cmd_str[cnt]):
                    mp = float(cmd_str[cnt])
                    self.tree.SetValue(MERGE_PARAMETER, mp)
                    self.operations.append(['b', mp])
                    print("MERGE_PARAMETER is set to " + str(mp))
                    self.cmd_history = self.cmd_history + ',b,' + str(mp)
            elif dat == 'v':
                self.tree.Validate()
            elif dat == 't':
                print(self.tree)
            elif dat == 'p':
                self.WriteToPython()
                self.cmd_history = self.cmd_history + ',p'
            elif dat == 'm':
                self.WriteToMatlab()
                self.cmd_history = self.cmd_history + ',m'
            elif dat == 'lv' or dat == 'lf' or dat == 'le': # list vertex or face or edge
                cnt = cnt + 1
                num = []
                while cnt < len(cmd_str):
                    if IsInt(cmd_str[cnt]):
                        num.append(int(cmd_str[cnt]))
                        cnt = cnt + 1
                    else:
                        cnt = cnt - 1
                        break
                if dat == "lf":
                    if num != []:
                        print("Faces to be listed: " + str(num))
                        for f in num:
                            self.tree.ListFace(f)
                    else:
                        self.tree.ListFaces()
                elif dat == "lv":
                    if num != []:
                        print("Vertices to be listed: " + str(num))
                        for v in num:
                            self.tree.ListVertex(v)
                    else:
                        self.tree.ListVertices()
                elif dat == "le":
                    self.tree.ListEdges()
            cnt = cnt + 1
        return 0

    def GetFileName(self):
        fn = self.name
        for op in self.operations:
            if op[0] == 'r':
                fn = fn + '_r'
                for f in op[1]:
                    fn = fn + '_' + str(f)
            elif op[0] == 'c':
                fn = fn + '_c'
                for f in op[1]:
                    fn = fn + '_' + str(f)
            elif op[0] == 'b':
                fn = fn + '_b'
                fn = fn + '_' + str(op[1])
        return fn.replace('.', 'p')

    def WriteToPython(self):
        fn = self.GetFileName()
        print("Print to python file " + fn + ".py")
        ifile = open(fn + ".py", 'w')
        ifile.write("import sys\n")
        ifile.write("import os\n")
        ifile.write("import " + self.name + "_include\n")
        ifile.write("import " + self.name + "_include\n")
        ifile.write("from " + self.name + "_include import *\n")
        ifile.write("model = " + self.name + "_include.Model('" + self.name + "',os.getcwd()+'/')\n")
        ifile.write("model.InitializeModel()\n")
        ifile.write("#command: " + self.cmd_history + "\n\n")
        ifile.write("tree = PolyTree2D()\n")
        ifile.write("tree.Synchronize(model.model_part, True)\n\n")
        for op in self.operations:
            if op[0] == 'r':
                for f in op[1]:
                    ifile.write("tree.MarkFaceRefine(" + str(f) + ")\n")
            elif op[0] == 'c':
                for f in op[1]:
                    ifile.write("tree.MarkFaceCoarsen(" + str(f) + ")\n")
            elif op[0] == 'b':
                ifile.write("tree.SetValue(MERGE_PARAMETER, " + str(op[1]) + ")\n\n")
            ifile.write("tree.BeginRefineCoarsen()\n")
            ifile.write("tree.EndRefineCoarsen()\n\n")
        ifile.write("tree.Validate()\n")
        ifile.write('tree.WriteMatlab("' + fn + '.m", True, True)\n')
        ifile.close()

    def WriteToMatlab(self):
        fn = self.GetFileName()
        print("Print to matlab file " + fn + ".m")
        self.tree.WriteMatlab(fn + ".m", True, True)
