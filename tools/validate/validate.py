import sys, msys
import unittest as UT

_mol = None

class TestCompatibility(UT.TestCase):
    def setUp(self): 
        self.mol=_mol

    def testContiguousGids(self):
        ''' Ensure that gids are contiguous, starting at zero.  '''
        gids=[a.gid for a in self.mol.atoms]
        self.assertEqual(gids, range(self.mol.natoms))

def Validate(mol):
    global _mol
    _mol=mol
    prg=UT.main(exit=False, argv=[sys.argv[0],'validate'])
    return prg.result


