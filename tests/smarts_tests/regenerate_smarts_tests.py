#!/usr/bin/garden-exec
#{
# garden env-keep-only --user
# garden load desres-python/3.7.2-08c7/bin
# garden prepend-path PYTHONPATH $(dirname $0)/../../build/lib/python
# python "$0" "$@"
#}

# Run from current directory (...msys/tests/smarts_tests)
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
        'PC2777373.mae',
        'PC32699.mae']

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
