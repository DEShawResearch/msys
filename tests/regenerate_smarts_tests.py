#!/usr/bin/env python2.7

import sys
import os
TMPDIR=os.getenv('TMPDIR', 'objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys
import json

files = ['/proj/desres/root/Linux/x86_64/dms_inputs/1.5.4/share/ww.dms',
        '/proj/desres/root/Linux/x86_64/dms_inputs/1.5.4/share/membrane.dms',
        'tests/files/4A9C_chainA_46.mae',
        'tests/files/4A9C_chainA_cf_21.mae',
        'tests/files/acrd.mae',
        'tests/files/azul.mae',
        'tests/files/boron.mae',
        'tests/files/fenz.mae',
        'tests/files/indz.mae',
        'tests/files/ndph.mae',
        'tests/files/pegm.mae']

mols = [msys.Load(f, structure_only=True) for f in files]
annot_mols = []
atoms = []
for mol in mols:
    msys.AssignBondOrderAndFormalCharge(mol)
    annot_mols.append(msys.AnnotatedSystem(mol))
    atoms.append(mol.select('not water'))

f = open('tests/smarts_tests_new.json', 'r')
tests = json.loads(f.read())
f.close()
for tup in tests:
    sp = msys.SmartsPattern(tup[0])
    for i in range(1, len(tup)):
        tup[i] = sp.findMatches(annot_mols[i-1], atoms[i-1])

json.dump(tests, file('tests/smarts_tests.json', 'w'))
