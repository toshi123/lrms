#!/bin/env python
# vim: set fileencoding=utf-8 :
#
# Author:   toshi123
# License:  MIT License
# Created:  2014-04-23
#

import sys
import argparse
import prody
from TMalign import TMalign


if __name__ == '__main__':
    usage = "usage: %prog reference_structure.pdb model_structure.pdb" 
    p = argparse.ArgumentParser(description="L-RMS calculater")
    p.add_argument('reference_PDBfile')
    p.add_argument('model_PDBfile')
    p.add_argument('-r','--ref_receptor',default='A',help='chain name of reference receptor')
    p.add_argument('-l','--ref_ligand',default='B',help='chain name of reference ligand')
    p.add_argument('-R','--model_receptor',default='A',help='chain name of model receptor')
    p.add_argument('-L','--model_ligand',default='B',help='chain name of model ligand')
    p.add_argument('--tmalign',default='/usr/local/bin/TMalign',help='path to TMalign')
    args = p.parse_args()

    ref_receptor = prody.parsePDB(args.reference_PDBfile,chain=args.ref_receptor)
    ref_ligand = prody.parsePDB(args.reference_PDBfile,chain=args.ref_ligand)
    model_receptor = prody.parsePDB(args.model_PDBfile,chain=args.model_receptor)
    model_ligand = prody.parsePDB(args.model_PDBfile,chain=args.model_ligand)

    tmalign = TMalign(model_receptor,ref_receptor,path = args.tmalign)
    trans = prody.Transformation(tmalign.matrix,tmalign.vector)
    trans.apply(model_ligand)
    lrms = prody.calcRMSD(model_ligand,ref_ligand)

    print lrms


