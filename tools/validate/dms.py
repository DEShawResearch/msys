from __future__ import print_function
from collections import defaultdict
import sys, msys
import unittest as UT
from msys import knot
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

    def testHasConstrainedHydrogens(self):
        ''' The hydrogens in the system have to be constrained '''
        constrained_hydrogens = []
        hydrogens = self.mol.select('hydrogen')
        hydrogenIds = [atom.id for atom in hydrogens]
        for table in self.mol.tables:
            if table.category=='constraint':
                constrained_hydrogens += [atom.id for term in table.findWithAny(hydrogens) for atom in term.atoms if atom.atomic_number == 1]
        unconstrained_hydrogens = list(set(hydrogenIds) - set(constrained_hydrogens))

        num = len(unconstrained_hydrogens)
        formatted_results = " ".join([str(aid) for aid in unconstrained_hydrogens[:10]]) + " ..." if len(unconstrained_hydrogens) > 10 else ""
        self.assertTrue(len(unconstrained_hydrogens) == 0, "Found %d hydrogens that are unconstrained: \n\t[%s]" % (num,
                formatted_results))

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
        ids = self.mol.selectIds('all')
        contacts = self.mol.findContactIds(1.0, ids, ignore_excluded=True)
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
            fragids = set([x.fragid for x in res.atoms])
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
        knots = knot.FindKnots(self.mol, max_cycle_size=10, ignore_excluded_knots=True)
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

    def testVirtualsInMultipleVirtualTables(self):
        ''' Virtuals should be in only one 'virtual_*' table '''
        # Ticket number Anton2#5176
        atomid_to_vtables = defaultdict(list)
        vtables = [t for t in self.mol.tables if t.name.startswith('virtual_')]
        vatoms = self.mol.select('atomicnumber 0')
        for table in vtables:
            for va in vatoms:
                terms = table.findWithAny([va])
                if len(terms) > 0:
                    atomid_to_vtables[va.id] += [table.name]
        virtuals_in_multiple_vtables = [aid for (aid, tables) in atomid_to_vtables.items()
                                       if len(tables) > 1]
        self.assertTrue(len(virtuals_in_multiple_vtables)==0, "Found virtuals that belong to multiple virtual_* tables: %s" % virtuals_in_multiple_vtables)


def Validate(mol, strict=False, verbose=1, anton=False,
        all=False):
    global _mol
    _mol=mol

    if isinstance(verbose, int):
        verbosity = verbose
    else:
        verbosity = 2 if verbose else 1

    if all:
        strict = True
    if strict:
        anton = True


    tests=[UT.TestLoader().loadTestsFromTestCase(TestBasic)]
    if anton:
        tests.append(UT.TestLoader().loadTestsFromTestCase(TestAnton))

    if strict: 
        tests.append(UT.TestLoader().loadTestsFromTestCase(TestStrict))

    suite=UT.TestSuite(tests)
    result=UT.TextTestRunner(verbosity=verbosity).run(suite)
    return result

main_doc = '''
dms-validate system.dms [ options ]

`dms-validate` flags conditions that are likely to be errors in a chemical
system.  The set of "basic" checks are always performed; additional checks
can be enabled using various command line flags. 

The set of basic checks comprise the following:

 * nonbonded: if a nonbonded table exists, every particle must have a 
   nonbonded param assignment.

 * knot: there must not be bonds passing through rings of size 10 or fewer.

 * nonzerobox: the volume of the global cell must be positive.
 
 * masslessatoms: there should be no massless atoms with atomic number > 0.

The set of strict checks comprise the following items.  Note that it
is certainly possible for a valid simulation to be performed using a
system that passes none of its strict checks!  However, it may be worth
investigating why a system fails theses checks.

 * constraints: the system must have constraint terms.  

 * consistent masses: Particles with equal atomic number must have equal mass.
   Pseudo particles (those with atomic number zero) are excluded from the
   check.

 * sparsify: every 1-4 bond (i.e., pair of atoms separated by three 
   distinct bonds) must be included in the exclusion table.

 * split waters: each water molecule must have its own unique resid,
   and live in its own residue.

'''

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ifile',
            help="Input file or jobid")
    parser.add_argument('--strict', action='store_true',
            help="Run strict tests")
    parser.add_argument('--verbose', action='store_true',
            help="Be chatty")
    parser.add_argument('--all', action='store_true',
            help="Run all tests")
    parser.add_argument('--anton', action='store_true',
            help="Run Anton-specific tests")
    return parser.parse_args()

def main():
    args = parse_args()
    mol = msys.Load(args.ifile)
    opts = args.__dict__
    opts.pop('ifile')
    result=Validate(mol, **opts)
    return not result.wasSuccessful()

