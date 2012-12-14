
import msys, sys, math

def count_overrides(nb):
    ''' return the number of pairwise nonbonded interactions which 
    are affected by the overrides in the given term table. '''
    # count the number of atoms for each nonbonded param
    nparams = nb.params.nparams
    pcount = [0] * nparams
    for t in nb.terms:
        pcount[t.param.id] += 1
    cnt = sum(pcount[pi.id] * pcount[pj.id] for pi, pj in nb.overrides())
    return cnt * 2
        
def duplicate_overrides(nb, pdict):
    # duplicate the overridden parameters and reassign while disambiguating
    nb.params.addProp('override', int)
    override_base = max(p['override'] for p in nb.params.params)+1
    for p, atoms in pdict.items():
        p_ = p.duplicate()
        p['override']=p.id+override_base
        p_['override']=p_.id+override_base
        for a in atoms: nb.term(a.id).param = p_
        # for the parameters that we duplicated, if there was already an
        # override involving that type, duplicate it with the new type.
        for (pi, pj), op in nb.overrides().items():
            if pi==p: pi=p_
            if pj==p: pj=p_
            nb.setOverride(pi,pj,op)

def add_override(nb, **params):
    p = nb.override_params.addParam()
    for k,v in params.items(): 
        nb.override_params.addProp(k,type(v))
        p[k]=v
    return p

def apply_override(nb, sel1, sel2, p):
    t1 = set(nb.term(a.id).param for a in sel1)
    t2 = set(nb.term(a.id).param for a in sel2)
    for pi in t1:
        for pj in t2:
            nb.setOverride(pi,pj,p)

def Override(mol, sel1, sel2, **params):
    ''' override the nonbonded interaction for selections sel1 and sel2
    with the given keyword parameters.  '''

    # add properties
    nb = mol.table('nonbonded')
    pdict = dict()
    for a in sel1+sel2:
        p = nb.term(a.id).param
        pdict.setdefault(p,[]).append(a)

    duplicate_overrides(nb, pdict)
    p = add_override(nb, **params)
    apply_override(nb, sel1, sel2, p)


def combine_arithmetic_geometric(pi,pj):
    sigma =        0.5*(pi['sigma'] + pj['sigma'])
    epsilon = math.sqrt(pi['epsilon'] * pj['epsilon'])
    return sigma, epsilon

def combine_geometric(pi,pj):
    sigma = math.sqrt(pi['sigma'] * pj['sigma'])
    epsilon = math.sqrt(pi['epsilon'] * pj['epsilon'])
    return sigma, epsilon

def Scale(mol, sel1, sel2, scale_sigma, scale_epsilon):
    ''' scale the sigma and epsilon interaction for selections s1 and s2
    with the given keyword parameters.  '''

    rule = mol.nonbonded_info.vdw_rule.lower()
    if rule=='arithmetic/geometric': 
        combine = combine_arithmetic_geometric
    elif rule=='geometric': 
        combine = combine_geometric
    else:
        raise ValueError, "Don't understand vdw_rule '%s'" % rule

    # partition selections by type, and call Override for each pair
    nb = mol.table('nonbonded')

    pdict1 = dict()
    pdict2 = dict()
    for a in sel1:
        p = nb.term(a.id).param
        pdict1.setdefault(p,[]).append(a)
    for a in sel2:
        p = nb.term(a.id).param
        pdict2.setdefault(p,[]).append(a)

    duplicate_overrides(nb, pdict1)
    duplicate_overrides(nb, pdict2)

    for p1, s1 in pdict1.items():
        for p2, s2 in pdict2.items():
            sigma, epsilon = combine(p1,p2)
            p = add_override(nb, 
                             sigma=scale_sigma*sigma, 
                             epsilon=scale_epsilon*epsilon)
            apply_override(nb, s1, s2, p)

