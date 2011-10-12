#!/usr/bin/env python2.7

import sys, unittest
sys.path.insert(0,'objs/Linux/x86_64/lib/python')
import msys

class TestMain(unittest.TestCase):

    def testAtom(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        self.assertTrue(a1.system==m)
        self.assertTrue(a2.system==m)
        self.assertEqual(a2.system.atoms[0], a1)
        self.assertEqual(a1.system.atoms[1], a2)
        self.assertNotEqual(m.atom(0), m.atom(1))

        # FIXME: hmm, maybe fetching m.atom(0) ought to throw after the
        # atom has been destroyed.  
        m.atom(0).destroy()
        self.assertFalse(m.atom(0) in m.atoms)

        m2=msys.CreateSystem()
        a3=m2.addAtom()
        self.assertNotEqual(a1,a3)
        self.assertEqual(a3,a3)
        self.assertEqual(a3,m2.atom(0))

        res=a1.residue
        self.assertEqual(m.residues[0], res)
        rnew = m.addResidue()
        self.assertNotEqual(rnew, a1.residue)
        a1.residue = rnew
        self.assertEqual(rnew, a1.residue)

    def testAtomProps(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        F, I, S = 'F', 'I', 'S'
        for p,t in zip((F, I, S), (float, int, str)):
            m.addAtomProp(p,t)
            self.assertTrue(p in m.atom_props)
            self.assertEqual(m.atomPropType(p), t)
        self.assertEqual(len(m.atom_props), 3)
        self.assertEqual(m.atom_prop_types, [float, int, str])
        a2=m.addAtom()

        a1.charge = 3.5
        self.assertEqual(a1.charge, 3.5)
        self.assertEqual(m.atom(a1.id).charge, a1.charge)
        self.assertEqual(a1[F], 0)
        self.assertEqual(a1[I], 0)
        self.assertEqual(a1[S], "")

        a1[F]=32.5
        a1[I]=42
        a1[S]="justinrocks"
        self.assertEqual(a1[F], 32.5)
        self.assertEqual(a1[I], 42)
        self.assertEqual(a1[S], "justinrocks")

        with self.assertRaises(KeyError):
            a1['foobar'] = 32

        self.assertTrue(hasattr(a2, 'charge'))
        self.assertTrue('F'in a2)
        self.assertFalse('foobar' in a2)

        m.delAtomProp('F')
        self.assertFalse('F'in a2)
        self.assertEqual(m.atom_prop_types, [int, str])


    def testBond(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        b=a2.addBond(a3)
        self.assertEqual(b.system, m)
        self.assertEqual(b.first, a2)
        self.assertEqual(b.second, a3)
        b.order=32.5
        self.assertEqual(b.order, 32.5)
        self.assertEqual(len(m.bonds),1)

        first, second = b.atoms
        self.assertEqual(first, a2)
        self.assertEqual(second, a3)

        self.assertEqual([a for a in a1.bonded_atoms], [])
        self.assertEqual([a for a in a2.bonded_atoms], [a3])
        self.assertEqual([a for a in a3.bonded_atoms], [a2])

        b.destroy()
        self.assertEqual(len(m.bonds),0)

        b1=a1.addBond(a2)
        self.assertEqual(len(a2.bonds),1)
        self.assertEqual(b1, m.findBond(a1,a2))
        self.assertEqual(b1, a1.findBond(a2))
        self.assertEqual(b1, a2.findBond(a1))
        self.assertEqual(None, m.findBond(a1,a3))

        b.destroy()
        b.destroy()
        b.destroy()
        self.assertEqual(len(a2.bonds),1)

        m.delBonds(m.bonds)
        self.assertEqual(len(m.bonds), 0)
        m.delAtoms(m.atoms)
        self.assertEqual(len(m.atoms), 0)

        self.assertEqual(len(m.residues), 3)
        m.delResidues(m.residues[2:])
        self.assertEqual(len(m.residues), 2)
        m.delResidues(m.residues[0:1])
        self.assertEqual(len(m.residues), 1)
        m.delChains(m.chains)
        self.assertEqual(len(m.chains), 0)

    def testRemove(self):
        m=msys.CreateSystem()
        r=m.addResidue()
        for i in range(10000):
            r.addAtom()
        r.destroy()
        r.destroy()
        r.destroy()
        self.assertEqual(len(m.atoms), 0)

        c=m.addChain()
        for i in range(10000):
            c.addResidue()
        c.destroy()
        c.destroy()
        c.destroy()
        self.assertEqual(len(m.residues), 0)

    def testResidue(self):
        m=msys.CreateSystem()
        r=m.addResidue()
        a1=r.addAtom()
        a2=r.addAtom()
        self.assertEqual(r.system, a1.system)
        self.assertEqual(r.system, m)
        self.assertEqual(m, a2.system)
        self.assertEqual(len(m.residues), 1)
        self.assertEqual(len(m.atoms), 2)
        self.assertEqual(len(r.atoms), 2)
        r.destroy()
        r.destroy()
        r.destroy()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(len(m.atoms), 0)
        self.assertEqual(len(r.atoms), 0)

    def testChain(self):
        m=msys.CreateSystem()
        c=m.addChain()
        r1=c.addResidue()
        r2=c.addResidue()
        a1=r1.addAtom()
        a2=r2.addAtom()
        self.assertEqual(c.system, r1.system)
        self.assertEqual(c.system, m)
        self.assertEqual(m, r2.system)
        self.assertEqual(len(m.chains), 1)
        self.assertEqual(len(m.residues), 2)
        self.assertEqual(len(c.residues), 2)
        self.assertEqual(len(m.atoms), 2)
        c.destroy()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(len(c.residues), 0)
        self.assertEqual(len(m.atoms), 0)

    def testTermTable(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        a4=m.addAtom()
        angle=m.addTable("angle", 3)
        self.assertEqual(angle.natoms, 3)
        self.assertEqual(len(angle.terms), 0)
        t1=angle.addTerm(m.atoms[:3])
        self.assertEqual(t1.atoms, m.atoms[:3])
        self.assertTrue(t1.param is None)
        t2=angle.addTerm(m.atoms[1:4], None)
        self.assertTrue(t2.param is None)
        self.assertEqual(len(angle.terms), 2)

        params=angle.param_table
        self.assertEqual(params, angle.param_table)
        p0=params.addParam()
        p1=params.addParam()
        self.assertEqual(p1.table, params)
        self.assertEqual(p0, params[0])
        self.assertEqual(p1, params[1])
        self.assertEqual(len(params), 2)
        t1.param=p1
        self.assertEqual(t1.param, p1)
        t1.param=None
        self.assertEqual(t1.param, None)

        with self.assertRaises(RuntimeError):
            angle.addTerm((m.atom(2),))
        with self.assertRaises(RuntimeError):
            angle.addTerm((m.atom(1),m.atom(1), m.atom(5)))

    def testDelTerm(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        a4=m.addAtom()
        table=m.addTable("foo", 2)
        table.addTerm((a1,a2))
        table.addTerm((a1,a3))
        table.addTerm((a1,a4))
        table.addTerm((a2,a3))
        table.addTerm((a4,a2))
        table.addTerm((a4,a3))
        self.assertEqual(table.nterms, 6)
        table.term(2).destroy()     # removes term 2
        self.assertEqual(table.nterms, 5)

        table.delTermsWithAtom(a3)  # removes term 1,3,5
        self.assertEqual(table.nterms, 2)
        a4.destroy()
        self.assertEqual(table.nterms, 1)


    def testParamProps(self):
        params=msys.CreateParamTable()
        self.assertEqual(len(params.props), 0)
        for n,t in zip(('F', 'I', 'S'), (float, int, str)):
            params.addProp(n,t)
            self.assertTrue(n in params.props)
            self.assertEqual(params.propType(n), t)
        self.assertEqual(params.props, ['F', 'I', 'S'])
        self.assertEqual(params.prop_types, [float, int, str])

        params.delProp('I')
        self.assertEqual(params.props, ['F', 'S'])
        self.assertEqual(params.prop_types, [float, str])

        p1=params.addParam()
        p1['F']=1.5
        p2=p1.duplicate()
        self.assertEqual(p1['F'], p2['F'])
        p2['F']=2.5
        self.assertNotEqual(p1['F'], p2['F'])


    def testDirectParamProps(self):
        ''' just like testParamProps, but operate directly on the TermTable
        instead of its ParamTable '''
        m=msys.CreateSystem()
        table=m.addTable("table", 5)
        self.assertEqual(len(table.term_props), 0)
        for n,t in zip(('F', 'I', 'S'), (float, int, str)):
            table.addTermProp(n,t)
            self.assertTrue(n in table.term_props)
            self.assertEqual(table.termPropType(n), t)

    def testParamWithNoProps(self):
        m=msys.CreateSystem()
        table=m.addTable("foo", 1)
        table.category='bond'
        self.assertEqual(table.category, 'bond')
        msys.SaveDMS(m, 'foo.dms')

    def testSharedParams(self):
        m1=msys.CreateSystem()
        m2=msys.CreateSystem()
        m1.addAtom()
        m2.addAtom()
        m2.addAtom()
        params=msys.CreateParamTable()
        p1=params.addParam()
        p2=params.addParam()
        table1=m1.addTable("table", 1, params)
        table2=m2.addTable("table", 1, params)
        self.assertEqual(table1.param_table, table2.param_table)
        t1=table1.addTerm(m1.atoms, p2)
        t2=table2.addTerm(m2.atoms[1:], p2)
        self.assertEqual(t1.param, t2.param)
        self.assertEqual(t2.param, p2)

        params.addProp("fc", float)
        p1['fc']=32
        p2['fc']=42
        self.assertEqual(t1.param['fc'],42)
        self.assertEqual(t2.param['fc'],42)

    def testTermProps(self):
        m=msys.CreateSystem()
        m.addAtom()
        table=m.addTable("table", 1)
        F="F"
        I="I"
        table.addTermProp(F, float)
        t1=table.addTerm(m.atoms)
        table.addTermProp(I, int)
        t2=table.addTerm(m.atoms)
        self.assertEqual(t1[F], 0)
        self.assertEqual(t1[I], 0)
        self.assertEqual(t2[F], 0)
        self.assertEqual(t2[I], 0)
        t1[F]=3.5
        t2[I]=42
        self.assertEqual(t1[F], 3.5)
        self.assertEqual(t1[I], 0)
        self.assertEqual(t2[F], 0)
        self.assertEqual(t2[I], 42)
        self.assertEqual(table.term_props, [F, I])
        self.assertEqual(table.term_prop_types, [float, int])

        m2=msys.CreateSystem()
        m2.addAtom()
        with self.assertRaises(RuntimeError):
            table.addTerm(m2.atoms)
        alist=m.atoms
        alist[0]=m2.atoms[0]
        with self.assertRaises(RuntimeError):
            table.addTerm(alist)
        alist[0]=m.atoms[0]
        table.addTerm(alist)

        table.delTermProp(I)
        with self.assertRaises(KeyError):
            t1[I]



    def testTableDMS(self):
        m=msys.CreateSystem()
        t=m.addTableFromSchema('stretch_harm')
        self.assertEqual(t.natoms, 2)
        self.assertEqual(t.props, ['r0', 'fc'])
        self.assertEqual(t.prop_types, [float, float])
        self.assertEqual(t.term_props, ['constrained'])
        self.assertEqual(t.termPropType('constrained'), int)

    def testNonbondedDMS(self):
        m=msys.CreateSystem()
        nb=m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        nb2=m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        self.assertEqual(nb,nb2)
        self.assertEqual(nb.props, ['sigma', 'epsilon'])
        self.assertEqual(nb.prop_types, [float,float])

    def testGlobalCell(self):
        m=msys.CreateSystem()
        m.cell.A.x=32
        self.assertEqual(m.cell.A.x, 32)
        tgt=[1,2,3]
        m.cell.A[:]=tgt
        m.cell[1][:]=tgt
        self.assertEqual(m.cell.A, m.cell.B)
        self.assertNotEqual(m.cell.A, m.cell.C)
        with self.assertRaises(IndexError): m.cell[3]
        with self.assertRaises(IndexError): m.cell[-4]

    def testClone(self):
        m=msys.CreateSystem()
        m.addAtom().name='a'
        m.addAtom().name='b'
        m.addAtom().name='c'
        b1=m.atom(0).addBond(m.atom(1))
        with self.assertRaises(KeyError):
            b1['res_order']=3.5
        m.addBondProp('res_order', float)
        b1['res_order']=3.5

        c=msys.CloneSystem(m.atoms[::-1])
        self.assertEqual( [a.name for a in c.atoms], ['c', 'b', 'a'])
        self.assertEqual(c.bond(0)['res_order'], 3.5)

        with self.assertRaises(RuntimeError):
            msys.CloneSystem([m.atoms[0], m.atoms[1], m.atoms[0]])

        m.atom(2).destroy()
        with self.assertRaises(RuntimeError):
            msys.CloneSystem([m.atom(0), m.atom(1), m.atom(2)])

    def testAppend(self):
        m=msys.CreateSystem()
        m.addAtom().name='a'
        m.addAtom().name='b'
        m.append(m)
        self.assertEqual( [a.name for a in m.atoms], ['a', 'b', 'a', 'b'])

    def testUpdateFragids(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        def fragids(m): return [a.fragid for a in m.atoms]

        frags=m.updateFragids()
        self.assertEqual(fragids(m), [0,1,2])
        self.assertEqual(frags, [[m.atom(i)] for i in range(3)])

        a1.addBond(a3)
        frags=m.updateFragids()
        self.assertEqual(fragids(m), [0,1,0])
        self.assertEqual(frags, [[m.atom(0), m.atom(2)], [m.atom(1)]])


if __name__=="__main__":
    unittest.main()
