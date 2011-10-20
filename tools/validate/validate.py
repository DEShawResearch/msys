import sys, msys
import unittest as UT

_mol = None

class TestCompatibility(UT.TestCase):
    def setUp(self): 
        self.mol=_mol

    def testContiguousGids(self):
        ''' Gids must be contiguous, starting at zero '''
        gids=[a.gid for a in self.mol.atoms]
        self.assertEqual(gids, range(self.mol.natoms))

class TestAnton(UT.TestCase):
    def setUp(self):
        self.mol=_mol

    def testHasConstraints(self):
        ''' The system must have constraints '''
        constraints=[t for t in self.mol.tables if t.category=='constraint']
        self.assertTrue(len(constraints)>0)

class TestSanity(UT.TestCase):
    def setUp(self):
        self.mol=_mol

    def testConsistentMasses(self):
        ''' Particles with equal atomic number must have equal mass '''
        m=dict()
        for a in self.mol.atoms:
            m.setdefault(a.atomic_number, set()).add(a.mass)
        for anum, masses in m.items():
            self.assertEqual(len(masses), 1, 
                    "Multiple masses for atomic number %s: \n\t%s" % (
                        anum, list(masses)))

def Validate(mol):
    global _mol
    _mol=mol
    prg=UT.main(exit=False, argv=[sys.argv[0],'validate'])
    return prg.result


