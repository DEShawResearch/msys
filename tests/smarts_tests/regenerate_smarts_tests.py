#!/usr/bin/env python2.7

# Run from current directory (...msys/tests/smarts_tests)
import sys
import os
TMPDIR=os.getenv('TMPDIR', '../../objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys
import ast
import pprint

files = ['ww.dms',
        'membrane.dms',
        '4A9C_chainA_46.mae',
        '4A9C_chainA_cf_21.mae',
        'acrd.mae',
        'azul.mae',
        'boron.mae',
        'fenz.mae',
        'indz.mae',
        'ndph.mae',
        'pegm.mae',
        'PC2777373.mae']

mols = []
annot_mols = []
atoms = []
names = []
for f in files:
    mols.append(msys.Load(f, structure_only=True))
    msys.AssignBondOrderAndFormalCharge(mols[-1])
    annot_mols.append(msys.AnnotatedSystem(mols[-1]))
    names.append(f[:-4])

for name, annot_mol in zip(names, annot_mols):
    f = open('%s_matches' % name, 'r')
    tests = ast.literal_eval(f.read())
    f.close()
    for k in tests.keys():
        sp = msys.SmartsPattern(k)
        tests[k] = sp.findMatches(annot_mol)
    pprint.pprint(tests, stream=open('%s_matches' % name, 'w'))
