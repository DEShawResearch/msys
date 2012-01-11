#!/usr/bin/env python2.7

import os, sys, unittest
TMPDIR=os.getenv('TMPDIR', 'objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys

from msys import builder

#defs=builder.Defs()
#defs.load('top_all27_prot_lipid_na.inp')
#
#mol=msys.LoadDMS('in.dms')
#builder.build(mol, defs)
#msys.SaveDMS(mol, 'out.dms')


# very high level - just specify paths
#builder.topology('top_all27_prot_lipid_na.inp')
#builder.build('in.dms', 'out.dms')

# next level - objects.  
defs = builder.Defs()
defs.load('top_all27_prot_lipid_na.inp')

mol=msys.LoadDMS('dhfr-in.dms', structure_only=True)
mol=mol.clone("noh")

# rename the one chain to 'A'
mol.chain(0).name='A'

# put all the water in its own chain
wat=mol.addChain()
wat.name="WAT"
for a in mol.select('water and noh'):
    a.residue.chain=wat

#mol=defs.build(mol) # uses default patches
defs.buildChain(mol.chain(0), plast='ACE')
defs.buildChain(wat)

defs.patch('HS2', mol.select('resid 44')[0].residue)

msys.SaveDMS(mol.sorted(), 'out.dms')

# next level - specific handling for chains.  Replace defs.build(mol) with

#defs.buildChain(mol.chain(0), first='ACE', last='CT3')
#defs.buildChain(mol.chain(1)) # water
#mol=mol.sorted()


