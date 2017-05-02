import sys, msys
import unittest as UT
from msys import knot, groups
from msys.wrap import Wrapper
import numpy

_mol = None

class TestCase(UT.TestCase):
    def setUp(self):
        self.mol=_mol

class TestAnton(TestCase):
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
                    "Use dms-fix-mass to fix multiple masses for atomic number %s: \n\t%s" % (
                        anum, list(masses)))

    def testOrthorhombicCell(self):
        ''' Unit cell must be diagonal '''
        cell=self.mol.getCell()
        self.assertTrue((numpy.diag(cell.diagonal())==cell).all(),
                "Unit cell is not diagonal")

    def testContacts(self):
        ''' No nonbonded atoms within 1A of each other '''
        ids = self.mol.selectIds('atomicnumber > 0')
        contacts = self.mol.findContactIds(1.0, ids)
        formatted_results = '\n'.join('%6d %6d %f' % x for x in contacts[:10])
        num = len(contacts)
        self.assertEqual(len(contacts), 0,
                "Found %d nonbonded atom pairs within 1A of each other: \n%s" % (num, formatted_results))

    def testPeriodicContacts(self):
        ''' No contacts crossing periodic boundaries within 1A '''
        mol=self.mol.clone()
        natoms=mol.natoms
        mol.append(mol)
        sel='(index < %d) and within 1 of (index >= %d)' % (natoms, natoms)
        pos = mol.positions
        a,b,c = mol.cell
        bad=[]
        for i in (-1,0,1):
            for j in (-1,0,1):
                for k in (-1,0,1):
                    if i==0 and i==j and i==k: continue
                    delta = i*a + j*b + k*c
                    pos[natoms:] += delta
                    mol.setPositions(pos)
                    ids = mol.selectIds(sel)
                    bad.extend(ids)
                    pos[natoms:] -= delta
        self.assertEqual(len(bad), 0,
                "Found atoms with periodic contacts less than 1A: %s" % str(bad))

class TestStrict(TestCase):
    def testSplitWaterResids(self):
        '''every water molecule must have its own resid'''
        water_atoms = self.mol.select('water and atomicnumber 8')

        help_string="Use dms-fix-water-residues to put each water in its own residue."

        for water_atom in water_atoms:
            res = water_atom.residue
            heavy = [a for a in res.atoms if a.atomic_number > 0]
            self.assertEqual(len(heavy), 3,
                             help_string)
            fragids = set(map(lambda x: x.fragid, res.atoms))
            self.assertEqual(len(fragids), 1,
                             help_string)

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
    def testKnots(self):
        ''' the system must not contain knots for rings of size <= 10 '''
        alchemical = self.mol.getTable('alchemical_nonbonded') is not None
        if alchemical:
            print "Skipping knot check on alchemical system"
        else:
            knots = knot.FindKnots(self.mol, max_cycle_size=10)
            self.assertTrue(len(knots)==0, "The system has %d bonds passing through small rings: %s" % (len(knots), knots))

    def testHasNonbonded(self):
        ''' Every particle must have a nonbonded param assignment '''
        self.assertTrue('nonbonded' in self.mol.table_names)
        nb=self.mol.table('nonbonded')
        for t in nb.terms:
            self.assertFalse(t.param is None,
                    "Particle with id %d has no nonbonded parameters" %
                    t.atoms[0].id)

    def testNonzeroBox(self):
        ''' The volume of the global cell must be positive. '''
        cell = self.mol.getCell()
        vol=numpy.linalg.det(cell)
        self.assertTrue(vol>0, "Global cell must have positive volume, got %s" % vol)

    def testMasslessAtoms(self):
        ''' There should be no massless atoms with atomic number > 0.'''
        massless_atoms = [atom.id for atom in self.mol.atoms if atom.atomic_number > 0 and atom.mass == 0]
        self.assertTrue(len(massless_atoms)==0, "Found massless atoms with atomic number > 0: %s" % massless_atoms)


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

    def testGroups(self):
        ''' Atoms in bonded terms must be within a clone radius of each
        other. '''
        radius = 5.5
        mol = self.mol.clone()
        Wrapper(mol).wrap()
        for table in mol.tables:
            self.assertFalse(groups.clone_buffer_violations(table, radius),
                    "Clone radius check failed for terms in '%s'; run dms-check-groups for more information." % table.name)
    

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

