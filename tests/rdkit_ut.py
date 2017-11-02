#!/usr/bin/garden-exec
#{
# source `dirname $0`/pyenv.sh
# garden load rdkit/2017.03.3-01c7/lib-python${PYSUFFIX}
# exec python $0 "$@"
#}
import os, sys, unittest
from util import *

import msys

import unittest
from rdkit import Chem
from time import time

class Main(unittest.TestCase):
    def testSmall(self):
        mol = msys.Load('tests/files/jandor-bad.sdf')
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

    def testChiralAtoms(self):
        mol = msys.Load('tests/files/jandor.sdf')
        rdmol = msys.ConvertToRdkit(mol)
        for r in rdmol.GetAtoms():
            if r.GetIdx() in (0,4,7):
                self.assertEqual(r.GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            elif r.GetIdx() in (5,):
                self.assertEqual(r.GetChiralTag(), Chem.ChiralType.CHI_TETRAHEDRAL_CCW)

    def testBondStereo(self):
        sdf = 'tests/files/34106.sdf'
        mol = msys.Load(sdf)
        rdmol = msys.ConvertToRdkit(mol)
        for i, r in enumerate(rdmol.GetBonds()):
            ai = r.GetBeginAtomIdx()
            aj = r.GetEndAtomIdx()
            if ai==16 and aj==17:
                self.assertEqual(r.GetStereo(), Chem.BondStereo.STEREOZ)
            else:
                self.assertEqual(r.GetStereo(), Chem.BondStereo.STEREONONE)

if __name__=="__main__":
  unittest.main(verbosity=2)
