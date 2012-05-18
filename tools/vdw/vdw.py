
import msys, sys

def Override(mol, s1, s2, **params):
    ''' override the nonbonded interaction for selections s1 and s2
    with the given keyword parameters.  '''

    # create combined table if needed
    cbname = 'nonbonded_combined'
    if not cbname in mol.table_names:
        cb = mol.addTable(cbname, 2)
        cb.category = 'override'
        for k,v in params.items():
            cb.params.addProp(k, type(v))
    else:
        cb = mol.table(cbname)
        assert cb.natoms == 2, "%s table must have 2 atoms per term, got %d" % (
                cbname, ct.natoms)

    sel1 = mol.select(s1)
    sel2 = mol.select(s2)
    print "Overriding interactions for %d x %d = %d particle pairs" % (
            len(sel1), len(sel2), len(sel1)*len(sel2))


    pdict = dict()
    nb = mol.table('nonbonded')
    for a in sel1 + sel2:
        p = nb.term(a.id).param
        pdict.setdefault(p,[]).append(a)

    # duplicate the overridden parameters and reassign while disambiguating
    nb.params.addProp('override', int)
    override_base = max(p['override'] for p in nb.params.params)+1
    for p, atoms in pdict.items():
        p_ = p.duplicate()
        p['override']=p.id+override_base
        p_['override']=p_.id+override_base
        for a in atoms: nb.term(a.id).param = p_

    # hash the existing tuples
    h=dict()
    for t in cb.terms:
        key = tuple(sorted(a.id for a in t.atoms))
        if key in h:
            raise ValueError, "duplicate tuple", t.atoms
        else:
            h[key] = t

    # add new tuples, replacing existing items
    p = cb.params.addParam()
    for k,v in params.items(): p[k] = v
    for a1 in sel1:
        for a2 in sel2:
            key = tuple(sorted((a1.id, a2.id)))
            term = h.get(key)
            if term is None:
                h[key] = cb.addTerm([a1,a2], p)
            else:
                term.param = p

