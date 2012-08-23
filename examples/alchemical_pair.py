#!/usr/bin/env python2.7

import msys
import os
import sys

ver=map(int, msys.version.split('.'))
assert ver >= (1,4,7), "Msys version 1.4.7 or later required."

def modify_dms(mol):

    if 'alchemical_pair_12_6_es' not in mol.table_names:
        print "WARNING: alchemical_pair_12_6_es table not found.  Nothing to do"
        return

    if 'alchemical_nonbonded' not in mol.table_names:
        print "WARNING: no alchemical particles.  Nothing to do."
        return
    
    # rename alchemical_pair_12_6_es to softcore.
    alcpair = mol.table('alchemical_pair_12_6_es')
    alcpair.name = 'alchemical_pair_softcore_es'

    # Add a chargeC property and initialize to zero.
    nb = mol.table('nonbonded')
    alc = mol.table('alchemical_nonbonded')
    alc.addTermProp('chargeC', float)
    for t in alc.terms: t['chargeC'] = 0.0

    # Alchemical atoms are 'core' if their nbtype doesn't change.
    core_terms = []
    for tB in alc.terms:
        a = tB.atoms[0]
        tA = nb.findExact([a])[0]
        paramA = tA.param
        paramB = tB.param
        # Set chargeC to average charge for core atoms.
        if paramA == paramB:
            tB['chargeC'] = 0.5*(tB['chargeB'] + a.charge)
            core_terms.append(tB)

    # select core-core pairs
    core_atoms = set(t.atoms[0] for t in core_terms)
    core_pairs = [t for t in alcpair.terms if set(t.atoms).issubset(core_atoms)]

    # Add core-core pairs to the regular pairs table
    pairs = mol.addTableFromSchema('pair_12_6_es')
    for t in core_pairs:
        p = pairs.params.addParam()
        qijA = t['qijA']
        qijB = t['qijB']
        qij  = 0.5*(qijA + qijB)
        p['aij'] = t['aijA']
        p['bij'] = t['bijA']
        p['qij'] = qij
        pairs.addTerm(t.atoms, p)

        # perform the alchemical charge operation if needed
        if qijA != qijB:
            t['aijA'] = 0.0
            t['aijB'] = 0.0
            t['bijA'] = 0.0
            t['bijB'] = 0.0
            t['qijA'] -= qij
            t['qijB'] -= qij
        else:
            t.remove()

if __name__=="__main__":
    if len(sys.argv) != 3:
        print "modify_dms.py input.dms output.dms"
        sys.exit(1)
    ipath= sys.argv[1]
    opath= sys.argv[2]

    try:
        os.unlink(opath)
    except OSError:
        pass
    mol=msys.Load(ipath)
    modify_dms(mol)
    mol.coalesceTables()
    mol=mol.clone()
    msys.SaveDMS(mol, opath)
