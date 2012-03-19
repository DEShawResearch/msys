#!/usr/bin/env desres-exec
#{
# desres-cleanenv -m Python/2.7.1-06A/bin -m msys/1.2.0/lib-python \
# -- python $0 "$@"
#}

import msys

def main(mol):
    # Add a chargeC property.  It will be nonzero only for alchemical atoms.
    mol.addAtomProp('chargeC', float)
    for a in mol.atoms:
        if a.alchemical:
            a['chargeC'] = 0.5*(a.charge + a.chargeB)

    # create a new pair_softcore_es table, with the same parameters as 
    # the pair_12_6_es table.
    hardpair = mol.table('pair_12_6_es')
    params = hardpair.params
    softpair = mol.addTable('pair_softcore_es', 2, params=params)
    softpair.category = 'bond'

    # modify the original alchemical hard pair terms
    for t in hardpair.terms:
        if not t.alchemical: continue
        A = t.param
        B = t.paramB
        # "noncore pair": if vdw parameters are mutating, move to softcore
        if A['aij']!=B['aij'] and A['bij']!=B['bij']:
            ts = softpair.addTerm(t.atoms)
            ts.param  = A
            ts.paramB = B
            t.remove()

        # "core pair": if vdw is unchanged but qij changes, create a new 
        # softcore pair term with some funny qij parameters and no vdw
        # parameters.  Modify the original pair term, and since all the
        # alchemical part is now handled by the new softcore term, remove
        # the alchemical part of the hardcore pair.
        elif A['aij']==B['aij'] and A['bij']==B['bij'] and A['qij']!=B['qij']:
            p  = A.duplicate()
            pA = params.addParam()
            pB = params.addParam()
            qC = (A['qij']+B['qij'])/2
            p ['qij'] = qC
            pA['qij'] = A['qij'] - qC
            pB['qij'] = B['qij'] - qC
            ts = softpair.addTerm(t.atoms)
            ts.param  = pA
            ts.paramB = pB
            t.param = p
            t.paramB = None

        # eliminate spurious alchemical pairs: those with identical A and B
        elif A['aij']==B['aij'] and A['bij']==B['bij'] and A['qij']==B['qij']:
            t.paramB = None

if __name__=="__main__":
    import sys
    ifile, ofile = sys.argv[1:]
    mol = msys.Load(ifile)
    main(mol)

    # the coalesce/clone trick gives us nonredundant parameter tables
    mol.coalesceTables()
    mol = mol.clone()
    msys.SaveDMS(mol, ofile)

