import sys, msys
import unittest as UT
from msys import knot


_mol = None

class TestCase(UT.TestCase):
    def setUp(self):
        self.mol=_mol

class TestAnton(TestCase):
    def testKnots(self):
        ''' the system must not contain knots for rings of size <= 10 '''
        knots = knot.FindKnots(self.mol, max_cycle_size=10)
        self.assertTrue(len(knots)==0, "The system has %d bonds passing through small rings: %s" % (len(knots), knots))

    def testHasConstraints(self):
        ''' The system must have constraints '''
        constraints=[t for t in self.mol.tables if t.category=='constraint']
        self.assertTrue(len(constraints)>0)

    def testConsistentMasses(self):
        ''' Particles with equal atomic number must have equal mass. 
        Pseudos excluded.  '''
        m=dict()
        for a in self.mol.atoms:
            m.setdefault(a.atomic_number, set()).add(a.mass)
        for anum, masses in m.items():
            if anum==0: continue
            self.assertEqual(len(masses), 1, 
                    "Multiple masses for atomic number %s: \n\t%s" % (
                        anum, list(masses)))

class TestStrict(TestCase):
    def testSparsify(self):
        ''' molecule must have all 1-4 exclusions.  '''
        # hash the exclusions and pairs
        excls=set()
        if not 'exclusion' in self.mol.table_names: return
        ptr=self.mol.table('exclusion')._ptr
        for t in ptr.terms():
            ai, aj = ptr.atoms(t)
            excls.add((ai,aj))
        p14=set()
        mol=self.mol._ptr
        for bi in mol.bonds():
            bnd=mol.bond(bi)
            b0 = bnd.i
            b1 = bnd.j
            for n0 in mol.bondedAtoms(b0):
                if n0==b1: continue
                for n1 in mol.bondedAtoms(b1):
                    if n1==b0: continue
                    if n1==n0: continue
                    pair=[n0, n1]
                    pair.sort()
                    p14.add(tuple(pair))

        ediff=p14.difference(excls)
        self.assertFalse(ediff, "%d 1-4 bonds not found in exclusions - is the file sparsified?" % len(ediff))

class TestBasic(TestCase):
    def testHasNonbonded(self):
        ''' Every particle must have a nonbonded param assignment '''
        self.assertTrue('nonbonded' in self.mol.table_names)
        nb=self.mol.table('nonbonded')
        for t in nb.terms:
            self.assertFalse(t.param is None,
                    "Particle with id %d has no nonbonded parameters" %
                    t.atoms[0].id)


class TestDesmond(TestCase):
    def testBondsBetweenNonbonded(self):
        ''' Flag bond terms or exclusions between atoms in different fragments.
        These atoms cannot be guaranteed to remain within a clone buffer
        radius of each other.  '''
        for table in self.mol.tables:
            if table.category in ('exclusion', 'bond') and table.natoms>1:
                for t in table.terms:
                    fragids = set(a.fragid for a in t.atoms)
                    self.assertEqual(len(fragids), 1, 
                            "Atoms %s form a term in the %s table but are not connected by bonds" % (t.atoms, table.name))
    

def Validate(mol, strict=False, desmond=False, verbose=False, anton=False,
        all=False):
    global _mol
    _mol=mol

    verbosity = 2 if verbose else 1
    if all:
        strict = True
        desmond = True
    if strict:
        anton = True


    tests=[UT.TestLoader().loadTestsFromTestCase(TestBasic)]
    if anton:
        tests.append(UT.TestLoader().loadTestsFromTestCase(TestAnton))

    if strict: 
        tests.append(UT.TestLoader().loadTestsFromTestCase(TestStrict))

    if desmond:
        tests.append(UT.TestLoader().loadTestsFromTestCase(TestDesmond))

    suite=UT.TestSuite(tests)
    result=UT.TextTestRunner(verbosity=verbosity).run(suite)
    return result

