#!/bin/env python
# vim: set fileencoding=utf-8 :
#
# Author:   toshi123
# License:  MIT License
# Created:  2014-03-25
#

import numpy as np
import tempfile
import subprocess
from prody import *

class TMalign:
    def __init__(self,pdbobj1,pdbobj2,path='/usr/local/bin/TMalign'):
        temp = tempfile.NamedTemporaryFile()
        pdbfile1 = self.writePDBtmp(pdbobj1)
        pdbfile2 = self.writePDBtmp(pdbobj2)
        tmalign_exe = [path,pdbfile1.name+".pdb",pdbfile2.name+".pdb","-m",temp.name]

        tmout = subprocess.check_output(tmalign_exe)
        matrixlines = temp.read().split("\n")[2:5]
        temp.close()
        pdbfile1.close()
        pdbfile2.close()

        lines = tmout.split("\n")
        self.score = self.getScore(lines)
        self.rmsd,self.seqid = self.getComparison(lines)
        self.vector,self.matrix = self.getRotationMatrix(matrixlines)

    def writePDBtmp(self,pdbobj):
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb")
        writePDB(tmp.name,pdbobj)
        return tmp

    def getScore(self,lines):
        length1 = int(lines[13].split()[3])
        length2 = int(lines[14].split()[3])
        score = 0
        if length1 < length2:
            score = float(lines[17].split()[1])
        else:
            score = float(lines[18].split()[1])
        return score

    def getComparison(self,lines):
        vals = lines[16].split(",")
        rmsd = float(vals[1].split()[1])
        seqid = float(vals[2].split()[1])
        return [rmsd,seqid]

    def getRotationMatrix(self,matrixlines):
        vec = []
        mat = []
        for line in matrixlines:
            vals = line.split()
            vec.append(vals[1])
            mat.append(vals[2:5])
        fvec = self.convertVectorString2Float(vec)
        fmat = self.convertMatrixString2Float(mat)
        vector = np.array(fvec)
        matrix = np.array(fmat) 
        return [vector,matrix]

    def convertMatrixString2Float(self,matrix):
        result = []
        for i in range(len(matrix)):
            result.append([])
            for j in range(len(matrix[i])):
                result[i].append(float(matrix[i][j]))
        return result

    def convertVectorString2Float(self,vector):
        result = []
        for i in range(len(vector)):
            result.append(float(vector[i]))
        return result


