#!/usr/bin/env python2.7

import os, sys, unittest
TMPDIR=os.getenv('TMPDIR', 'objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys
import numpy as NP

def vsize():
    cmd='ps -p %d h -o vsz' % os.getpid()
    s=os.popen(cmd).read()
    return int(s)


class TestMain(unittest.TestCase):

    def testSSSR(self):
        sys = msys.LoadDMS('/d/en/gregerse-0/p4/sw/forcefields/viparr4/cubane.dms', True)
        rings = msys.GetSSSR(sys.atoms, False)
        self.assertTrue(len(rings) == 5)
        for ring in rings:
            self.assertTrue(len(ring) == 4)
            for ring2 in rings:
                if ring != ring2:
                    intersect = set([a.id for a in ring]) & set([a.id for a in ring2])
                    self.assertTrue(len(intersect) == 2 or len(intersect) == 0)
        rings = msys.GetSSSR(sys.atoms, True)
        self.assertTrue(len(rings) == 6)
        for ring in rings:
            self.assertTrue(len(ring) == 4)
            for ring2 in rings:
                if ring != ring2:
                    intersect = set([a.id for a in ring]) & set([a.id for a in ring2])
                    self.assertTrue(len(intersect) == 2 or len(intersect) == 0)

    def testAssignBondOrdersAndFormalCharges(self):
        # Smoke test only
        sys = msys.LoadDMS('/proj/desres/root/Linux/x86_64/dms_inputs/1.5.4/share/ww.dms')
        msys.AssignBondOrderAndFormalCharge(sys)
        self.assertFalse('resonant_charge' in sys.atom_props)
        self.assertFalse('resonant_order' in sys.bond_props)
        msys.AssignBondOrderAndFormalCharge(sys.select('water'))
        self.assertFalse('resonant_charge' in sys.atom_props)
        self.assertFalse('resonant_order' in sys.bond_props)

    def testMemoryLeakPosVel(self):
        mol=msys.CreateSystem()
        for i in xrange(1000):
            mol.addAtom()

        pos=NP.zeros((mol.natoms, 3), 'f')
        for i in xrange(100):
            mol.setPositions(pos)
            mol.setVelocities(pos)
        oldsize = vsize()
        for i in xrange(1000):
            mol.setPositions(pos)
            mol.setVelocities(pos)
        newsize = vsize()
        self.assertEqual(oldsize, newsize)

    def testBadPath(self):
        m=msys.CreateSystem()
        m.addAtom()
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m,'/root/nono.dms')
        with self.assertRaises(RuntimeError):
            msys.SaveMAE(m,'/root/nono.mae')
        with self.assertRaises(RuntimeError):
            msys.SavePDB(m,'/root/nono.pdb')

    def xxtestFailedWrite(self):
        print "building big system..."
        m=msys.CreateSystem()
        c=m.addChain()
        r=c.addResidue()
        id=r.id
        f=r._ptr.addAtom
        for i in range(1000000):
            f(id)
        print "saving..."
        msys.SaveMAE(m,'/usr/tmp/big.dms')
        print "done"

    def testAtom(self):

        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        self.assertTrue(a1.system==m)
        self.assertTrue(a2.system==m)
        self.assertEqual(a2.system.atoms[0], a1)
        self.assertEqual(a1.system.atoms[1], a2)
        self.assertNotEqual(m.atom(0), m.atom(1))

        self.assertNotEqual(a1,None)
        self.assertTrue(a1!=None)
        self.assertTrue(a1==a1)
        self.assertFalse(a1==None)
        self.assertFalse(a1!=a1)


        # FIXME: hmm, maybe fetching m.atom(0) ought to throw after the
        # atom has been removed.  
        m.atom(0).remove()
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

    def testGlue(self):
        m=msys.CreateSystem()
        m.addAtom()
        m.addAtom()
        m.addAtom()

        with self.assertRaises(RuntimeError): m.addGluePair(1,1)
        with self.assertRaises(RuntimeError): m.addGluePair(1,10)
        with self.assertRaises(RuntimeError): m.addGluePair(10,1)

        self.assertEqual(m.gluePairs(), [])
        self.assertFalse(m.hasGluePair(1,2))
        self.assertFalse(m.hasGluePair(2,1))
        m.addGluePair(2,1)
        m.addGluePair(0,1)
        self.assertTrue(m.hasGluePair(1,0))
        self.assertTrue(m.hasGluePair(1,2))
        self.assertTrue(m.hasGluePair(2,1))
        self.assertEqual(m.gluePairs(), [(0,1),(1,2)])
        m2=m.clone('index 1 2')
        m3=m2.clone()
        m3.append(m)
        m.delGluePair(2,1)
        self.assertFalse(m.hasGluePair(1,2))
        self.assertFalse(m.hasGluePair(2,1))
        self.assertEqual(m2.gluePairs(), [(0,1)])
        self.assertEqual(m3.gluePairs(), [(0,1), (2,3),(3,4)])

    def testSelectOnRemoved(self):
        m=msys.CreateSystem()
        m.addAtom()
        m.addAtom().atomic_number=17
        m.addAtom()
        for a in m.select('atomicnumber 17'): a.remove()
        self.assertEqual(0, len(m.select('atomicnumber 17')))

    def testParamHandle(self):
        params=msys.CreateParamTable()
        self.assertNotEqual(params,None)
        self.assertEqual(params,params)
        self.assertNotEqual(params, msys.CreateParamTable())
        p1=params.addParam()
        p2=params.addParam()
        self.assertTrue(p1==p1)
        self.assertTrue(p1!=p2)
        self.assertFalse(p1!=p1)
        self.assertFalse(p1==p2)
        p1b = params.param(0)
        self.assertTrue(p1==p1b)
        self.assertTrue(p1b!=p2)
        self.assertFalse(p1b!=p1)
        self.assertFalse(p1b==p2)

        self.assertEqual(len(set((p1,p2))), 2)
        self.assertEqual(len(set((p1,p1b))), 1)

    def testSystemHandle(self):
        m=msys.CreateSystem()
        self.assertNotEqual(m,None)
        self.assertEqual(m,m)
        self.assertNotEqual(m,msys.CreateSystem())

    def testTermHandle(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        T=m.addTable('foo', 1)
        t1=T.addTerm([a], None)
        t2=T.addTerm([a], None)
        t1b=T.term(0)

        self.assertEqual(len(set((t1,t2))), 2)
        self.assertEqual(len(set((t1,t1b))), 1)

    def testCoalesce(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        table=m.addTableFromSchema('posre_harm')
        p1=table.params.addParam()
        p2=table.params.addParam()
        p3=table.params.addParam()
        p3['fcx'] = 32


        t1=table.addTerm([a], p1)
        t2=table.addTerm([a], p2)
        t3=table.addTerm([a], p3)

        self.assertFalse( t1.param==t2.param)
        self.assertFalse( t1.param==t3.param)

        table.coalesce()
        self.assertTrue( t1.param==t2.param)
        self.assertTrue( t1.param==p1)
        self.assertTrue( t2.param==p1)
        self.assertFalse( t1.param==t3.param)
        self.assertTrue( t3.param==p3)


    def testAtomProps(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        F, I, S = 'F', 'I', 'S'
        for p,t in zip((F, I, S), (float, int, str)):
            m.addAtomProp(p,t)
            self.assertTrue(p in m.atom_props)
            self.assertEqual(m.atomPropType(p), t)
        self.assertEqual(len(m.atom_props), 3)
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
        self.assertEqual([m.atomPropType(n) for n in m.atom_props], 
                [int, str])


    def testBond(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        b=a2.addBond(a3)

        self.assertTrue(a1!=b)
        self.assertFalse(a1==b)

        self.assertEqual(b.system, m)
        self.assertEqual(b.first, a2)
        self.assertEqual(b.second, a3)
        b.order=32
        b.resonant_order=32.5
        self.assertEqual(b.order, 32)
        self.assertEqual(b.resonant_order, 32.5)
        self.assertEqual(len(m.bonds),1)

        first, second = b.atoms
        self.assertEqual(first, a2)
        self.assertEqual(second, a3)

        self.assertEqual([a for a in a1.bonded_atoms], [])
        self.assertEqual([a for a in a2.bonded_atoms], [a3])
        self.assertEqual([a for a in a3.bonded_atoms], [a2])

        b.remove()
        self.assertEqual(len(m.bonds),0)

        b1=a1.addBond(a2)
        self.assertEqual(len(a2.bonds),1)
        self.assertEqual(b1, m.findBond(a1,a2))
        self.assertEqual(b1, a1.findBond(a2))
        self.assertEqual(b1, a2.findBond(a1))
        self.assertEqual(None, m.findBond(a1,a3))

        self.assertEqual(a2.nbonds, 1)

        b.remove()
        b.remove()
        b.remove()
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
        self.assertEqual(r.natoms, 10000)
        r.remove()
        r.remove()
        r.remove()
        self.assertEqual(len(m.atoms), 0)

        c=m.addChain()
        for i in range(10000):
            c.addResidue()
        self.assertEqual(c.nresidues, 10000)

        c.remove()
        c.remove()
        c.remove()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(c.nresidues, 0)

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
        r.name='alA'
        r.resid=32
        self.assertEqual(r.name, 'alA')
        self.assertEqual(r.resid, 32)
        
        r.remove()
        r.remove()
        r.remove()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(len(m.atoms), 0)
        self.assertEqual(len(r.atoms), 0)

    def testSegid(self):
        m=msys.CreateSystem()
        with self.assertRaises(RuntimeError):
            m.addAtomProp('segid', str)
        r=m.addResidue()
        a1=r.addAtom()
        a2=r.addAtom()
        r.chain.segid="WAT1"
        msys.SaveDMS(m,'foo.dms')
        m2=msys.LoadDMS('foo.dms')
        msys.Load("foo.dms")
        self.assertEqual(1,m.nresidues)
        self.assertEqual(1,m2.nresidues)

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
        c.remove()
        self.assertEqual(len(m.residues), 0)
        self.assertEqual(len(c.residues), 0)
        self.assertEqual(len(m.atoms), 0)

    def testResetParams(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        T=m.addTable('foo', 1)
        P=T.params
        P.addProp('f', float)
        p=P.addParam()
        p['f']=32
        t=T.addTerm([a1], p)

        P2=msys.CreateParamTable()
        P2.addProp('i', int)
        p2=P2.addParam()
        p2['i']=4

        with self.assertRaises(RuntimeError):
            T.term(0).param = p2

        T.params = P2
        self.assertEqual(T.params, P2)
        self.assertEqual(T.term(0).param, None)
        T.term(0).param = p2
        with self.assertRaises(RuntimeError):
            T.term(0).param = p

        T.params = P
        self.assertEqual(T.params, P)
        T.term(0).param = p


    def testTermTable(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        a4=m.addAtom()
        with self.assertRaises(RuntimeError):
            m.addTable("foo", 0)

        angle=m.addTable("angle", 3)
        self.assertEqual(angle.natoms, 3)
        self.assertEqual(len(angle.terms), 0)
        t1=angle.addTerm(m.atoms[:3])

        self.assertTrue(angle.hasTerm(0))
        self.assertFalse(angle.hasTerm(1))

        self.assertEqual(t1.atoms, m.atoms[:3])
        self.assertTrue(t1.param is None)
        t2=angle.addTerm(m.atoms[1:4], None)
        self.assertTrue(t2.param is None)
        self.assertEqual(len(angle.terms), 2)

        angle.name = 'angle2'
        self.assertEqual(angle.name, 'angle2')
        self.assertTrue('angle2' in m.table_names)
        self.assertFalse('angle' in m.table_names)
        angle.name = 'angle'
        self.assertEqual(angle.name, 'angle')
        self.assertTrue('angle' in m.table_names)
        self.assertFalse('angle2' in m.table_names)

        params=angle.params
        self.assertEqual(params, angle.params)
        p0=params.addParam()
        p1=params.addParam()
        self.assertEqual(p1.table, params)
        self.assertEqual(p0, params.param(0))
        self.assertEqual(p1, params.param(1))
        self.assertEqual(params.nparams, 2)
        t1.param=p1
        self.assertEqual(t1.param, p1)
        t1.param=None
        self.assertEqual(t1.param, None)

        with self.assertRaises(RuntimeError):
            angle.addTerm((m.atom(2),))
        with self.assertRaises(RuntimeError):
            angle.addTerm((m.atom(1),m.atom(1), m.atom(5)))

        angle.remove()
        self.assertFalse('angle' in m.table_names)
        self.assertFalse(angle in m.tables)

    def testIndex(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        a3=m.addAtom()
        a4=m.addAtom()
        t=m.addTable("foo", 2)

        self.assertEqual(t.findWithAll([a1]), [])
        self.assertEqual(t.findWithAny([a1]), [])

        t.addTerm((a1,a2))
        t.addTerm((a1,a3))
        t.addTerm((a1,a4))
        t.addTerm((a2,a3))
        t.addTerm((a4,a2))
        t.addTerm((a4,a3))
        self.assertEqual(t.nterms, 6)
        self.assertEqual(t.findWithAll([a1]), [t.term(x) for x in (0,1,2)])
        self.assertEqual(t.findWithAll([a2]), [t.term(x) for x in (0,3,4)])
        self.assertEqual(t.findWithAll([a3]), [t.term(x) for x in (1,3,5)])
        self.assertEqual(t.findWithAll([a4]), [t.term(x) for x in (2,4,5)])
        self.assertEqual(t.findWithAny([a1,a2]), [t.term(x) for x in range(5)])
        self.assertEqual(t.findWithAny([a1,a3]), [t.term(x) for x in (0,1,2,3,5)])
        self.assertEqual(t.findWithAll([a4,a2]), [t.term(4)])
        self.assertEqual(t.findWithAll([a2,a4]), [t.term(4)])
        self.assertEqual(t.findExact([a4,a2]), [t.term(4)])
        self.assertEqual(t.findExact([a2,a4]), [])
        a1.remove()
        self.assertEqual(t.findWithAll([a1]), [t.term(x) for x in ()])
        self.assertEqual(t.findWithAll([a2]), [t.term(x) for x in (3,4)])
        self.assertEqual(t.findWithAll([a3]), [t.term(x) for x in (3,5)])
        self.assertEqual(t.findWithAll([a4]), [t.term(x) for x in (4,5)])
        a5=m.addAtom()
        t.addTerm((a5,a3))
        self.assertEqual(t.findWithAll([a5]), [t.term(6)])
        self.assertEqual(t.findWithAll([a3]), [t.term(x) for x in (3,5,6)])
        self.assertEqual(t.findWithAny([a3,a5]), [t.term(x) for x in (3,5,6)])

    def testParamIndex(self):
        params=msys.CreateParamTable()
        params.addProp('i', int)
        params.addProp('f', float)
        params.addProp('s', str)
        p0=params.addParam()
        p1=params.addParam()
        p2=params.addParam()
        p3=params.addParam()
        p0['i']=1
        p1['i']=1
        p1['f']=1.5
        p2['f']=1.5
        p2['s']='x'
        p3['s']='x'
        for (k,v), ids in {
                ('i',0) : (2,3),
                ('i',1) : (0,1),
                ('f',0) : (0,3),
                ('f',1.5) : (1,2),
                ('s','') : (0,1),
                ('s','x') : (2,3)
                }.items():
            self.assertEqual(params.find(k,v), [params.param(x) for x in ids])
        for i in range(100): params.addParam()
        for (k,v), ids in {
                ('i',1) : (0,1),
                ('f',1.5) : (1,2),
                ('s','x') : (2,3)
                }.items():
            self.assertEqual(params.find(k,v), [params.param(x) for x in ids])
        p3['s']='y'
        for (k,v), ids in {
                ('i',1) : (0,1),
                ('f',1.5) : (1,2),
                ('s','x') : (2,)
                }.items():
            self.assertEqual(params.find(k,v), [params.param(x) for x in ids])

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
        table.term(2).remove()     # removes term 2
        self.assertEqual(table.nterms, 5)

        table.delTermsWithAtom(a3)  # removes term 1,3,5
        self.assertEqual(table.nterms, 2)
        a4.remove()
        self.assertEqual(table.nterms, 1)


    def testParamProps(self):
        params=msys.CreateParamTable()
        self.assertEqual(len(params.props), 0)
        for n,t in zip(('F', 'I', 'S'), (float, int, str)):
            params.addProp(n,t)
            self.assertTrue(n in params.props)
            self.assertEqual(params.propType(n), t)
        self.assertEqual(params.props, ['F', 'I', 'S'])

        params.delProp('I')
        self.assertEqual(params.props, ['F', 'S'])

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

    def testPropOverlap(self):
        ''' ensure term props cannot overlap param props. '''
        m=msys.CreateSystem()
        table=m.addTable("table", 1)
        table.params.addProp("a", float)
        with self.assertRaises(RuntimeError):
            table.addTermProp("a", int)


    def testParamWithNoProps(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        table=m.addTable("foo", 1)
        p=table.params.addParam()
        table.addTerm([a], p)
        table.category='bond'
        self.assertEqual(table.category, 'bond')
        msys.SaveDMS(m, 'foo.dms')

    def testFunnyNames(self):
        m=msys.CreateSystem()
        aux=msys.CreateParamTable()
        aux.addProp('group', float)
        aux.addParam()
        m.addAuxTable('aux', aux)
        msys.SaveDMS(m, 'foo.dms')

        m2=m.clone()
        aux2=m2.auxtable('aux')
        self.assertEqual(aux.nprops, aux2.nprops)
        self.assertEqual(aux.nparams, aux2.nparams)

        m3=msys.CreateSystem()
        m3.append(m)
        aux3=m3.auxtable('aux')
        self.assertEqual(aux.nprops, aux3.nprops)
        self.assertEqual(aux.nparams, aux3.nparams)


        table=m.addTable('foo', 1)
        table.category='bond'
        table.addTermProp('select', int)
        table.params.addProp('group', float)
        m.addAtom()
        table.addTerm(m.atoms, table.params.addParam())
        msys.SaveDMS(m, 'bar.dms')

    def testFormalCharge(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        self.assertEqual(a.formal_charge, 0)
        a.formal_charge=32
        self.assertEqual(a.formal_charge, 32)
        a.formal_charge=-10
        self.assertEqual(a.formal_charge, -10)
        msys.SaveDMS(m, 'bar.dms')
        m2=msys.LoadDMS('bar.dms')
        self.assertEqual(m.atom(0).formal_charge, -10)
                
    def testRefcount(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        table=m.addTable("foo", 1)
        ptr=table.params._ptr
        p1=table.params.addParam()
        p2=table.params.addParam()
        p3=table.params.addParam()
        t1=table.addTerm([a], p1)
        self.assertEqual(ptr.refcount(t1.param.id), 1)
        t2=table.addTerm([a], p1)
        self.assertEqual(ptr.refcount(t1.param.id), 2)
        self.assertEqual(ptr.refcount(t2.param.id), 2)
        t1.remove()
        self.assertEqual(ptr.refcount(t2.param.id), 1)
        t3=table.addTerm([a], p2)
        self.assertEqual(ptr.refcount(t2.param.id), 1)
        self.assertEqual(ptr.refcount(t3.param.id), 1)

        t3.remove()
        self.assertEqual(ptr.refcount(t2.param.id), 1)


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
        self.assertEqual(table1.params, table2.params)
        t1=table1.addTerm(m1.atoms, p2)
        t2=table2.addTerm(m2.atoms[1:], p2)
        self.assertEqual(t1.param, t2.param)
        self.assertEqual(t2.param, p2)

        params.addProp("fc", float)
        p1['fc']=32
        p2['fc']=42
        self.assertEqual(t1.param['fc'],42)
        self.assertEqual(t2.param['fc'],42)

        # refcounts work even across shared tables
        t1['fc']=52
        self.assertEqual(t2['fc'], 42)

    def testTermParamProps(self):
        m=msys.CreateSystem()
        m.addAtom()
        table=m.addTable("table", 1)
        p=table.params.addParam()
        t=table.addTerm(m.atoms, p)
        t2=table.addTerm(m.atoms, p)

        with self.assertRaises(KeyError):
            t['x']
            t.stateB['x']
            t.stateA['x']

        table.params.addProp('x', float)
        table.params.addProp('y', str)
        self.assertEqual(t['x'], 0)
        self.assertEqual(t['y'], '')

        t2['x']=32
        t2['y']='whodat'
        self.assertEqual(t['x'], 0)
        self.assertEqual(t['y'], '')
        self.assertEqual(table.params.nparams, 2)

        self.assertEqual(t2._ptr.getProp(t2.id, 0), 32)
        t2._ptr.setProp(t2.id, 0, 33)
        self.assertEqual(t2._ptr.getProp(t2.id, 0), 33)
        t2._ptr.setProp(t2.id, 0, 32)

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

    def testAux(self):
        m=msys.CreateSystem()
        e=msys.CreateParamTable()
        m.addAuxTable("cmap", e)
        self.assertTrue('cmap' in m.auxtable_names)
        self.assertEqual(e, m.auxtable('cmap'))
        self.assertTrue(e in m.auxtables)
        m.delAuxTable('cmap')
        self.assertFalse('cmap' in m.auxtable_names)
        self.assertFalse(e in m.auxtables)

    def testSchemas(self):
        m=msys.CreateSystem()
        for s in msys.TableSchemas():
            m.addTableFromSchema(s)
        for s in msys.NonbondedSchemas():
            m.addNonbondedFromSchema(s)
            m.nonbonded_info.vdw_funct=""

    def testTableDMS(self):
        assert 'stretch_harm' in msys.TableSchemas()
        assert 'vdw_12_6' in msys.NonbondedSchemas()

        m=msys.CreateSystem()
        t=m.addTableFromSchema('stretch_harm')
        self.assertEqual(t.natoms, 2)
        self.assertEqual(t.params.props, ['r0', 'fc'])
        self.assertEqual([t.params.propType(n) for n in t.params.props], 
                [float, float])
        self.assertEqual(t.term_props, ['constrained'])
        self.assertEqual(t.termPropType('constrained'), int)

    def testNonbondedDMS(self):
        m=msys.CreateSystem()
        nb=m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        nb2=m.addNonbondedFromSchema("vdw_12_6", "arithemetic/geometric")
        m.addNonbondedFromSchema("vdw_12_6")
        self.assertEqual(nb,nb2)
        params=nb.params
        props=params.props
        self.assertEqual(props, ['sigma', 'epsilon'])
        self.assertEqual([nb.params.propType(n) for n in props], [float,float])

    def testMultipleNonbonded(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        m.addAtom()
        disp = m.addTable('nonbonded_dispersion', 1)
        repl = m.addTable('nonbonded_repulsion', 1)
        elec = m.addTable('nonbonded_charge', 1)
        for t in disp, repl, elec: 
            t.category='nonbonded'
            p=t.params.addParam()
            t.addTerm([a],p)

        disp.params.addProp('foo', float)
        repl.params.addProp('bar', float)
        elec.params.addProp('charge', float)
        m.nonbonded_info.vdw_funct = "disp_repl_charge"
        m.nonbonded_info.vdw_rule = "geom/geom/geom"

        fname='/tmp/saveit.dms'
        try:
            msys.SaveDMS(m,fname)
        finally:
            pass
            #if os.path.exists(fname):
                #os.unlink(fname)
        t=m.addTable('nonbonded', 1)
        t.category='nonbonded'
        p=t.params.addParam()
        t.addTerm([a], p)
        t.addTerm([m.atom(1)], p)
        try:
            msys.SaveDMS(m,fname)
        finally:
            pass
        


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

        # accept numpy floats instead of doubles. 
        c=NP.zeros((3,3), 'f')
        m.setCell(c)
        # FIXME: this is hard to make work due to how I implemented the
        # bindings.  
        #m.cell.A[:] = c[0]

    def testNonbondedInfo(self):
        m=msys.CreateSystem()
        nb=m.nonbonded_info
        nb.vdw_funct = "justinrocks"
        nb.vdw_rule = "yep"
        self.assertEqual(nb.vdw_funct , "justinrocks")
        self.assertEqual(nb.vdw_rule , "yep")

    def testClone2(self):
        m=msys.CreateSystem()
        a1=m.addAtom()
        a2=m.addAtom()
        t1=m.addTable('t1', 1)
        t1.category='nonbonded'
        params=t1.params
        t2=m.addTable('t2', 1, params)
        t2.category='nonbonded'

        p1=params.addParam()
        p2=params.addParam()
        t1.addTerm([a1], p1)
        t1.addTerm([a2], p2)
        t2.addTerm([a2], p1)

        self.assertEqual(params.nparams, 2)

        m2=m.clone()
        t1=m2.table('t1')
        t2=m2.table('t2')
        params=t1.params
        self.assertEqual(t1.params,t2.params)
        self.assertEqual(params.nparams, 2)

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

        m.atom(2).remove()
        with self.assertRaises(RuntimeError):
            msys.CloneSystem([m.atom(0), m.atom(1), m.atom(2)])

    def testAppend(self):
        m=msys.CreateSystem()
        a0=m.addAtom()
        a0.name='a'
        a1=m.addAtom()
        a1.name='b'
        b=a0.addBond(a1)
        m.addAtomProp('foo', float)
        m.addBondProp('bar', int)
        a1['foo']=3.14
        b['bar']=42
        m.append(m)
        self.assertEqual( [a.name for a in m.atoms], ['a', 'b', 'a', 'b'])
        self.assertEqual( m.atom(2)['foo'], 0 )
        self.assertEqual( m.atom(3)['foo'], 3.14 )
        self.assertEqual( m.bond(1)['bar'], 42)

    def testRemoveTable(self):
        m=msys.CreateSystem()
        a=m.addAtom()
        T=m.addTable('foo', 1)
        T.params.addProp('x', int)
        p=T.params.addParam()
        t1=T.addTerm([a], p)
        t1['x']=32
        self.assertEqual(t1.param.id,0)

        T2=m.addTable('bar', 1, T.params)
        t2=T2.addTerm([a], p)
        t2['x']=32
        self.assertEqual(t2.param.id,1)
        self.assertEqual(T.params.nparams,2)

        T.name= 'fooj'
        self.assertEqual(T.system, m)
        self.assertEqual(T.name, 'fooj')
        T.remove()
        self.assertEqual(T.system, None)
        self.assertEqual(T.name, "")
        with self.assertRaises(RuntimeError):
            T.addTerm([a], p)
        self.assertEqual(T.name, '')
        self.assertEqual([t.id for t in T.terms], [])

    def testChargeC(self):
        m=msys.CreateSystem()
        nb=m.addNonbondedFromSchema('vdw_12_6')
        alc=m.addTable('alchemical_nonbonded', 1, nb.params)
        alc.category='nonbonded'
        alc.addTermProp('chargeC', float)

        a=m.addAtom()
        p=nb.params.addParam()

        nb.addTerm([a], p)
        t=alc.addTerm([a], p)
        t['chargeC']=0.5

        msys.SaveDMS(m, 'foo.dms')
        m2=msys.LoadDMS('foo.dms')
        msys.SaveDMS(m2, 'foo.dms')
        alc=m2.table('alchemical_nonbonded')
        self.assertEqual(alc.term(0)['chargeC'], 0.5)


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

    def testBadDMS(self):
        path='/tmp/_tmp_.dms'
        m=msys.CreateSystem()
        msys.SaveDMS(m, path)
        a0=m.addAtom()
        msys.SaveDMS(m, path)
        nb=m.addNonbondedFromSchema('vdw_12_6', 'geometric')
        # can't save - no terms for particle 0
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        t=nb.addTerm([a0])
        # can't save - no param
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        p0=nb.params.addParam()
        t.param=p0
        # now we can save
        msys.SaveDMS(m, path)
        # twiddle param to something bad
        t._ptr.setParam(t.id, None)
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        t._ptr.setParam(t.id, 0)

        stretch=m.addTableFromSchema('stretch_harm')
        t=stretch.addTerm([a0,a0])
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)
        p=stretch.params.addParam()
        t.param=p
        msys.SaveDMS(m, path)
        t._ptr.setParam(t.id, None)
        with self.assertRaises(RuntimeError):
            msys.SaveDMS(m, path)

    def testPositions(self):
        m=msys.CreateSystem()
        a0=m.addAtom()
        a1=m.addAtom()
        a2=m.addAtom()
        a0.y=1
        a1.z=2
        a2.x=3
        p=m.positions
        self.assertEqual(p[0][1],1)
        self.assertEqual(p[1][2],2)
        self.assertEqual(p[2][0],3)
        p[1][2]=4
        m.positions=p
        p=m.positions
        self.assertEqual(p[1][2],4)
        self.assertTrue((p==m.positions).all())
        with self.assertRaises(ValueError):
            m.positions=p[1:,:]
        with self.assertRaises(ValueError):
            m.positions=p[:,1:]

        # test of the low level Id-based accessors for positions
        ids=m._ptr.select('index 1')
        m._ptr.setPositions([[8,9,10]], ids)
        self.assertEqual(list(m.positions[1]), [8,9,10])
        self.assertEqual(list(m._ptr.getPositions(ids)[0]), [8,9,10])

    def testVelocities(self):
        m=msys.CreateSystem()
        a0=m.addAtom()
        a1=m.addAtom()
        a2=m.addAtom()
        a0.vy=1
        a1.vz=2
        a2.vx=3
        p=m.getVelocities()
        self.assertEqual(p[0][1],1)
        self.assertEqual(p[1][2],2)
        self.assertEqual(p[2][0],3)
        p[1][2]=4
        m.setVelocities(p)
        p=m.getVelocities()
        self.assertEqual(p[1][2],4)
        self.assertTrue((p==m.getVelocities()).all())
        with self.assertRaises(ValueError):
            m.setVelocities(p[1:,:])
        with self.assertRaises(ValueError):
            m.setVelocities(p[:,1:])

        # test of the low level Id-based accessors for velocities
        ids=m._ptr.select('index 1')
        m._ptr.setVelocities([[8,9,10]], ids)
        self.assertEqual(list(m.getVelocities()[1]), [8,9,10])
        self.assertEqual(list(m._ptr.getVelocities(ids)[0]), [8,9,10])

    def testMacros(self):
        m=msys.CreateSystem()
        m.addAtom().name='CA'
        m.addAtom().name='CB'
        self.assertFalse('foobar' in m.selection_macros)
        with self.assertRaises(RuntimeError): m.select('foobar')
        m.addSelectionMacro('foobar', 'name CB')
        self.assertEqual(m.selectionMacroDefinition('foobar'), 'name CB')
        self.assertTrue('foobar' in m.selection_macros)
        self.assertEqual(m.select('foobar')[0].id, 1)

        m2=m.clone()
        self.assertEqual(m2.select('foobar')[0].id, 1)

        m.delSelectionMacro('foobar')
        with self.assertRaises(RuntimeError): m.select('foobar')

        m.addSelectionMacro('foo', 'name CA')
        m.addSelectionMacro('bar', 'foo')
        m.addSelectionMacro('foo', 'bar')
        with self.assertRaises(RuntimeError):
            m.select('foo')

    def testAlchemicalMaeRestraint(self):
        d=os.path.dirname(__file__)
        m=msys.Load(os.path.join(d, 'mae/alchemical_restraint.mae'))
        self.assertEqual(m.natoms, 5)
        con=m.table('constraint_ah2')
        res=m.table('posre_harm')
        self.assertEqual(con.nterms, 1)
        self.assertEqual(sorted(x.id for x in con.term(0).atoms), [0,1,2])
        self.assertEqual(res.nterms, 2)
        self.assertEqual(res.term(0).atoms[0].id, 0)
        self.assertEqual(res.term(1).atoms[0].id, 1)
        self.assertEqual(res.term(0)['fcx'], 0.25)
        self.assertEqual(res.term(1)['fcy'], 0.35)

if __name__=="__main__":
    unittest.main()
