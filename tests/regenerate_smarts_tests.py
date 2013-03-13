#!/usr/bin/env python2.7

import sys
import os
TMPDIR=os.getenv('TMPDIR', 'objs/Linux/x86_64')
sys.path.insert(0,os.path.join(TMPDIR, 'lib', 'python'))
import msys
import json

mols = [msys.Load('/proj/desres/root/Linux/x86_64/dms_inputs/1.5.4/share/ww.dms'), msys.Load('/proj/desres/root/Linux/x86_64/dms_inputs/1.5.4/share/membrane.dms')]
annot_mols = []
atoms = []
matches = []
for mol in mols:
    msys.AssignBondOrderAndFormalCharge(mol)
    annot_mols.append(msys.AnnotatedSystem(mol))
    atoms.append(mol.select('not water'))
    matches.append([])

f = open('tests/smarts_tests.json', 'r')
d = json.loads(f.read())
f.close()
for patt in d['smarts']:
    sp = msys.SmartsPattern(patt)
    for i in range(2):
        matches[i].append(sp.findMatches(annot_mols[i], atoms[i]))

d['ww_matches'] = matches[0]
d['membrane_matches'] = matches[1]
json.dump(d, file('tests/smarts_tests.json', 'w'), separators=(', \n', ': '))
