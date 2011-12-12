
def apply(mol, atoms, fcx, fcy, fcz, replace=False):
    ''' add position restraints to atoms '''

    table=mol.addTableFromSchema('posre_harm')
    index={}
    if replace:
        for t in table.terms: t.remove()
    else:
        for t in table.terms:
            index[t.atoms[0]] = t

    param=table.params.addParam()
    param['fcx']=fcx
    param['fcy']=fcy
    param['fcz']=fcz

    for a in atoms:
        t=index.get(a)
        if t is None:
            t=table.addTerm([a])
        t.param = param
        t['x0']=a.x
        t['y0']=a.y
        t['z0']=a.z

    return table.nterms
