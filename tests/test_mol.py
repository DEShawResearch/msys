#!/usr/bin/env python2.7

import sys
sys.path.insert(0,'objs/lib/python')

import msys

ent=msys.CreateSystem()

atm=ent.addAtom()
print "len(atm.bonds):", len(atm.bonds)

print "old pos:", atm.pos
atm.pos=[1,2,3]
atm.pos=(1,2,.3)
print "new pos:", atm.pos

a2=ent.addAtom()
b=atm.addBond(a2)
print atm, a2, b.first, b.second

for bond in atm.bonds:
    print "bond for atm %s: %s" % (atm, bond)
del bond

ai=b.first
aj=b.second
print ai, aj, ai.id, aj.id

assert b.first == b.first
assert b.first != b.second

print "b id:", b.id
print "atm id:", atm.id
assert b != atm

atm3=ent.addAtom()
print "atm3 id:", atm3.id
assert not atm3 == atm

print "residue:", atm.residue
assert len(atm.residue.atoms)==1
assert len(atm.residue.chain.residues)==1

print "ent:", ent

for a in ent.atoms:
    print "ent %s atom %s" % (ent, a.id)

print [a for a in ent.atoms]
print "I have %d atoms" % len(ent.atoms)
print "Here's the last atom:", ent.atoms[-1]



del ai
del aj
del atm
del a2
del ent
del a

import gc
gc.collect()

print b
print b.first
print b.first
print b.first

del b

gc.collect()
print "collected!"

print dir()

for path in sys.argv[1:]:
    c=msys.ImportDMS(path)
    print "%d atoms, %d bonds, %d residues, %d chains" % (
            len(c.atoms), len(c.bonds), len(c.residues), len(c.chains))
    for name in c.table_names:
        table=c.table(name)
        natoms=table.natoms
        terms=table.terms
        p=table.params
        t=terms[0]
        print "  %20s: %8d terms,  %d atoms" % (name,len(terms),natoms)
        print "  %20s: %8d params, %d props" % ("",   len(p), len(p.props))
        print "  first term: ", t.id, t.param, [x.id for x in t.atoms]
        p0 = t.param
        for n in p.props: print "      %20s: %s" % (n, p0[n])
        if 'ff' in p.props:
            p0['ff'] = 32
        if 'fc' in p.props:
            p0['fc'] = 3.14
        for n in p.props: print "      %20s: %s" % (n, p0[n])

        # switch up a term
        if 'fc' in p.props and len(table.params)>1:
            print "old fc:", t.param['fc']
            t.param = table.params[1]
            print "new fc:", t.param['fc']

        assert len([t for t in terms])==len(terms)


