#!/bin/env python
# vim: set fileencoding=utf-8 :
#
# Author:   toshi123
# License:  MIT License
# Created:  2014-04-23
#

import sys
from optparse import OptionParser
import prody

if __name__ == '__main__':
    usage = "usage: %prog reference_structure.pdb model_structure.pdb" 
    p = OptionParser(usage=usage)
    (options, args ) = p.parse_args()
    if len(args) != 2:
        p.error( "incorrect number of arguments" )
    reffile,modelfile = args

    ref = prody.parsePDB(reffile)
    model = prody.parsePDB(modelfile)

    length_a = len(ref.select('chain A and ca').getResnums())
    length_b = len(ref.select('chain B and ca').getResnums())
    if length_a > length_b:
        receptor = 'chain A'
        ligand = 'chain B'
    else:
        receptor = 'chain B'
        ligand = 'chain A'

    r_matches = prody.matchChains(ref.select(receptor),model.select(receptor), pwalign=True)[0]
    l_matches = prody.matchChains(ref.select(ligand),model.select(ligand), pwalign=True)[0]

    t = prody.calcTransformation(r_matches[0], r_matches[1])
    t.apply(model)
    lrmsd = prody.calcRMSD(l_matches[0],l_matches[1])

    print lrmsd


