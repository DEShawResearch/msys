
import msys, sys

def Override(mol, s1, s2, **params):
    ''' override the nonbonded interaction for selections s1 and s2
    with the given keyword parameters.  '''

    # add properties
    nb = mol.table('nonbonded')
    cb = nb.overrides

    # perform selections
    sel1 = mol.select(s1)
    sel2 = mol.select(s2)

    pdict = dict()
    for a in sel1 + sel2:
        p = nb.term(a.id).param
        pdict.setdefault(p,[]).append(a)

    # duplicate the overridden parameters and reassign while disambiguating
    nb.params.addProp('override', int)
    override_base = max(p['override'] for p in nb.params.params)+1
    L = cb.list()
    for p, atoms in pdict.items():
        p_ = p.duplicate()
        p['override']=p.id+override_base
        p_['override']=p_.id+override_base
        for a in atoms: nb.term(a.id).param = p_
        # for the parameters that we duplicated, if there was already an
        # override involving that type, duplicate it with the new type.
        for i,j in L:
            if i==j and p.id==i:
                key=(p_.id, p_.id)
            elif p.id==i:
                key=(p_.id, j)
            elif p.id==j:
                key=(i, p_.id)
            else:
                continue
            cb.set(key, cb.get((i,j)))

    p = nb.override_params.addParam()
    for k,v in params.items(): 
        nb.override_params.addProp(k,type(v))
        p[k]=v

    # get the new set of types
    t1 = set(nb.term(a.id).param.id for a in sel1)
    t2 = set(nb.term(a.id).param.id for a in sel2)
    for i in t1:
        for j in t2:
            cb.set((i,j),p.id)
    
