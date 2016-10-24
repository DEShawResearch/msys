#!/usr/bin/garden-exec
#{
# source `dirname $0`/../MODULES
# exec garden with -c -e TMPDIR -m $PYTHON/bin \
# -m RDKit/2016_03_1-04/lib-python \
# -- python $0 "$@"
#}

import os, sys
TMPDIR=os.getenv('TMPDIR', 'objs/%s/x86_64' % os.getenv('DESRES_OS'))
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys

import unittest
from rdkit import Chem
from time import time

class Main(unittest.TestCase):
    def testSmall(self):
        mol = msys.Load('tests/files/jandor.sdf')
        with self.assertRaises(ValueError):
            rdmol = msys.ConvertToRdkit(mol)
        msys.AssignBondOrderAndFormalCharge(mol)
        amol = msys.AnnotatedSystem(mol)
        rdmol = msys.ConvertToRdkit(mol)
        for matm, ratm in zip(mol.atoms, rdmol.GetAtoms()):
            self.assertEqual(matm.atomic_number, ratm.GetAtomicNum())
            self.assertEqual(matm.formal_charge, ratm.GetFormalCharge())
            self.assertEqual(amol.valence(matm), ratm.GetExplicitValence())

    def testBig(self):
        mol = msys.Load('tests/files/2f4k.dms')
        msys.AssignBondOrderAndFormalCharge(mol)
        t=-time()
        rdmol = msys.ConvertToRdkit(mol)
        t+=time()
        print "%s: %d atoms, %d bonds in %.3fs" % (mol.name, mol.natoms, mol.nbonds, t)



if __name__=="__main__":
  unittest.main(verbosity=2)
