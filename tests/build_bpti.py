#!/usr/bin/env python2.7

import os, sys, unittest
TMPDIR=os.getenv('TMPDIR', 'objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys

from msys import builder

defs = builder.Defs()
defs.load('top_all27_prot_lipid_na.inp')

# 5pti saved as pdb from VMD
mol=msys.LoadPDB('5pti.pdb')

# ignore everything except protein and crystal water
mol=mol.clone('resid 1 to 58 or resname DOD')

# put the water in its own chain
wat=mol.addChain()
wat.name='XWAT'
for a in mol.select('resname DOD'): a.residue.chain=wat

# rename residues to match topology definitions
for old, new in {
        'CYS' : 'CYSH',
        'LYS' : 'LYSH',
        'DOD' : 'SOL',
        }.items():
    for a in mol.select('resname %s' % old):
        a.residue.name=new

# rename atom names to match topology definitions
for (oldres, old), new in {
        ('SOL', 'O')    : 'OH2',
        ('SOL', 'D1')   : 'H1',
        ('SOL', 'D2')   : 'H2',
        }.items():
    sel=mol.select('resname %s and name %s' % (oldres, old))
    print "renaming %d atoms from %s:%s to %s" % (len(sel), oldres, old, new)
    for a in sel: a.name=new

# build it
for c in mol.chains: defs.buildChain(c)

# add disulfide bridges
defs.patch('DISU', *[a.residue for a in mol.select('name CA and resid 14 38')])
defs.patch('DISU', *[a.residue for a in mol.select('name CA and resid 5 55')])
defs.patch('DISU', *[a.residue for a in mol.select('name CA and resid 30 51')])

# center it
mol.translate([-p for p in mol.center])

# solvate it
from msys import solvate
mol=solvate.Solvate(mol, dims=[64])

msys.SaveDMS(mol.sorted(), 'build.dms')


