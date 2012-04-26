
import _msys
from msys import kept
from math import sqrt

def find_residue(mol, resname, resid, chainname):
    for chn in mol.chains():
        if mol.chain(chn).name == chainname:
            for res in mol.residuesForChain(chn):
                residue = mol.residue(res)
                if residue.name == resname and residue.resid == resid:
                    return chn, res
            return chn, -1
    return -1, -1

def canonical(ids):
    n=len(ids)
    if   n==2: a,b=ids
    elif n==3: a,b=ids[0],ids[2]
    elif n==4: a,b=ids[1],ids[2]
    else: raise ValueError, "Expect len(ids) in (2,3,4), got %d" % len(ids)
    if a>b:
        ids.reverse()
    return ids

def copy_table( mol, btable ):
    atable = mol.addTable(btable.name(), btable.atomCount())
    ap=atable.params()
    bp=btable.params()
    for i in range(bp.propCount()):
        ap.addProp(bp.propName(i), bp.propType(i))
    for i in range(bp.termPropCount()):
        ap.addTermProp(bp.termPropName(i), bp.termPropType(i))
    return atable


def make_block( mol, amap, bmap, apairs, bpairs, atable, btable):
    if atable is None and btable is None:
        return []

    if atable is None:
        atable = copy_table(mol, btable)

    # a_item_list: the ids of the atoms in each atom, canonicalized
    natoms = atable.atomCount()
    arange = range(natoms)
    block=[]

    # bterms: mapping from canonicalized atoms to 1-based index
    bterms={}
    for i in range(btable.termCount()) if btable is not None else []:
        ids=btable.atoms(i)
        key=canonical([bmap[j] for j in ids])
        bterms[tuple(key)]=i+1

    for i in range(atable.termCount()):
        ids=atable.atoms(i)
        key=canonical([amap[j] for j in ids])
        j = bterms.pop(tuple(key), 0)
        block.append( [[i+1, j], key] )

    # add remaining items in b_item_dict
    block.extend( sorted( [[0,j], list(ids)] for ids, j in bterms.items() ))

    # Assign kept (-1) and not kept (0) for terms involving only dummy atoms.
    for inds, term in block:
        ti, tj = inds
        if tj<1:
            dums=[x for x in term if bpairs[x]<0]
            inds[1] = -1 if len(dums)==len(term) else 0
        elif ti<1:
            dums=[x for x in term if apairs[x]<0]
            inds[0] = -1 if len(dums)==len(term) else 0

    return block

def make_pairmaps(mol, bmap, atable, btable):
    if atable is None and btable is None:
        return []

    if atable is None:
        atable = copy_table(mol, btable)

    # make a and b items same way as exclmap
    b_item_dict=dict()
    for t in btable.terms():
        item = sorted(bmap[a] for a in btable.atoms(t))
        b_item_dict[tuple(item)]=t+1

    block=list()
    for t in atable.terms():
        item = sorted(a for a in atable.atoms(t))
        j=b_item_dict.get(tuple(item), 0) # 0 instead of -1 as in exclusion
        block.append( [[t+1,j], item] )
        if j!=0:
            del b_item_dict[tuple(item)]

    block2=list()
    for item, j in b_item_dict.items():
        block2.append( [[0,j], item] )
    block2.sort()
    block.extend(block2)

    return block

def get_pairs_table(mol):
    funct = mol.nonbonded_info.vdw_funct
    if funct=='vdw_12_6':
        pname = 'pair_12_6_es'
    elif funct=='vdw_exp_6':
        pname = 'pair_exp_6_es'
    else:
        raise ValueError, "Unsupported vdw_funct '%s'" % funct
    return mol.table(pname)

def find_nonbonded_param(table, atom):
    for t in table.terms():
        if table.atom(t,0)==atom:
            return table.param(t)
    raise ValueError, "Missing nonbonded term for atom %s" % atom

def convert_sig_eps(sij, eij):
    aij = pow(sij,12) * eij * 4.0
    bij = pow(sij, 6) * eij * 4.0
    return aij, bij

def combine_geometric(vi, vj):
    sij = sqrt(vi[0] * vj[0])
    eij = sqrt(vi[1] * vj[1])
    return convert_sig_eps(sij, eij)

def combine_arith_geom(vi,vj):
    sij = 0.5*(vi[0] + vj[0])
    eij = sqrt(vi[1] * vj[1])
    return convert_sig_eps(sij, eij)

def configure_vdw_12_6_param(pairs, param, nb, rule, nb_i, nb_j):
    ''' Using the given nonbonded entries, configure the vdw part of
    the given pair entry. '''
    # get sigma, epsilon for each atom
    isig = nb.propIndex('sigma')
    ieps = nb.propIndex('epsilon')
    iprops = [nb.getProp(nb_i, isig), nb.getProp(nb_i, ieps)]
    jprops = [nb.getProp(nb_j, isig), nb.getProp(nb_j, ieps)]

    if rule=='geometric':
        aij, bij = combine_geometric(iprops, jprops)
    elif rule=='arithmetic/geometric':
        aij, bij = combine_arith_geom(iprops, jprops)
    else:
        raise ValueError, "unsupported rule '%s' for vdw_12_6" % rule
    iaij = pairs.propIndex('aij')
    ibij = pairs.propIndex('bij')
    pairs.setProp(param, iaij, aij)
    pairs.setProp(param, ibij, bij)

def make_full_pair_entry(mol, ai, aj):
    pairs = get_pairs_table(mol)
    params = pairs.params()
    param = params.addParam()
    index = params.propIndex('qij')
    params.setProp(param, index, mol.atom(ai).charge * mol.atom(aj).charge)

    nb = mol.table('nonbonded')
    nb_i = find_nonbonded_param(nb, ai)
    nb_j = find_nonbonded_param(nb, aj)
    funct = mol.nonbonded_info.vdw_funct
    rule = mol.nonbonded_info.vdw_rule

    if funct=='vdw_12_6':
        configure_vdw_12_6_param( params, param, nb.params(), rule, nb_i, nb_j )
    elif funct=='vdw_exp_6':
        configure_vdw_exp_6_param( params, param, nb.params(), rule, nb_i, nb_j )
    else:
        raise ValueError, "Unsupported vdw_funct '%s'" % funct

    return pairs, param

def make_exclmap(bmap, apairs, bpairs, atable, btable):
    invbmap=dict()
    for key, val in bmap.items(): invbmap[val]=key

    pairs = get_pairs_table(atable.system())
    name = "alchemical_"+pairs.name()
    alcpairs = atable.system().addTableFromSchema(name,name)

    b_item_dict=dict()
    for t in btable.terms():
        item = sorted(bmap[a] for a in btable.atoms(t))
        b_item_dict[tuple(item)]=t+1

    block=list()
    for t in atable.terms():
        item=sorted(a for a in atable.atoms(t))
        j=b_item_dict.get(tuple(item), -1)
        block.append( [[t+1,j], item] )
        if j!=-1:
            del b_item_dict[tuple(item)]
            continue
        # An excluded pair of atoms in A mapped onto an unexcluded pair of
        # atoms in B.  We keep the exclusion, but add an extra pair to B.
        ai=invbmap.get(item[0], -1)
        aj=invbmap.get(item[1], -1)
        if apairs[item[0]]>=0 and apairs[item[1]]>=0 and ai>=0 and aj>=0:
            #print 'Offset exclusion %d-%d in A by pair interaction in B: %d-%d' % (apairs[item[0]], apairs[item[1]], ai, aj)

            pairsB, paramB = make_full_pair_entry(btable.system(), ai, aj)
            # remove non-alchemical version if it exists
            for t in pairs.terms():
                if sorted(pairs.atoms(t))==item:
                    pairs.delTerm(t)
                    break

            # find or create the alchemical pair
            for t in alcpairs.terms():
                if sorted(alcpairs.atoms(t))==item:
                    p = alcpairs.param(t)
                    break
            else:
                ids=_msys.IdList()
                for i in item: ids.append(i)
                p = alcpairs.params().addParam()
                t = alcpairs.addTerm(ids, p)

            # assign its alchemical state
            for i in range(pairsB.params().propCount()):
                prop = pairsB.params().propName(i)
                col = alcpairs.params().propIndex(prop+'B')
                alcpairs.params().setProp(p,col,pairsB.params().getProp(paramB,i))

    block2=list()
    for item, j in b_item_dict.items():
        block2.append( [[-1,j], item] )
        ai=invbmap.get(item[0], -1)
        aj=invbmap.get(item[1], -1)
        if apairs[item[0]]>=0 and apairs[item[1]]>=0 and ai>=0 and aj>=0:
            #print 'Offset exclusion %d-%d in B by pair interaction in A: %d-%d' % (ai, aj, apairs[item[0]], apairs[item[1]])
            pairs,param = make_full_pair_entry(atable.system(),item[0],item[1])
            # remove non-alchemical version if it exists
            for t in pairs.terms():
                if sorted(pairs.atoms(t))==sorted(item):
                    pairs.delTerm(t)
                    break

            # find or create the alchemical pair
            for t in alcpairs.terms():
                if sorted(alcpairs.atoms(t))==sorted(item):
                    p = alcpairs.param(t)
                    break
            else:
                ids=_msys.IdList()
                for i in item: ids.append(i)
                p = alcpairs.params().addParam()
                t = alcpairs.addTerm(ids, p)

            # assign its alchemical state
            for i in range(pairs.params().propCount()):
                prop = pairs.params().propName(i)
                col = alcpairs.params().propIndex(prop+'A')
                alcpairs.params().setProp(p,col,pairs.params().getProp(param,i))

    block2.sort()
    block.extend(block2)
    return block

def add_more_exclusions( pairs, bmap, A, B, cutoff=6.0 ):
    r2=cutoff*cutoff
    a_list=list()
    b_list=list()
    block=list()
    xyz=('x', 'y', 'z')

    for i, (a,b) in enumerate(pairs):
        if a>=0 and b<0: a_list.append(a)
        if a<0 and b>=0: b_list.append(b)

    for a in a_list:
        atm=A.atom(a)
        apos=[getattr(atm, x) for x in xyz]
        for b in b_list:
            atm=B.atom(b)
            bpos=[getattr(atm, x) for x in xyz]
            d2=sum((bpos[i] - apos[i])**2 for i in range(3))
            if d2<r2:
                block.append( [[-1,-1], [a,bmap[b]]] )
                #print "extra:", a, bmap[b]

    return block


def copy_param( dstparams, srcparams, srcid ):
    ''' copy parameters from param with id srcid in srctable to dsttable '''
    if dstparams == srcparams:
        return srcparams.duplicate(srcid)
    dstid = dstparams.addParam()
    for dstindex in range(dstparams.propCount()):
        prop = dstparams.propName(dstindex)
        srcindex = srcparams.propIndex(prop)
        if srcindex<0: continue
        newval = srcparams.getProp(srcid, srcindex)
        dstparams.setProp(dstid, dstindex, newval)
    return dstid

def copy_alchemical(dst, srcA, idA, srcB, idB):
    ''' copy params from src into alchemical dst '''
    id = dst.addParam()
    if not _msys.bad(idA):
        for i in range(srcA.propCount()):
            prop = srcA.propName(i)
            col = dst.propIndex(prop+'A')
            dst.setProp(id,col,srcA.getProp(idA,i))
    if not _msys.bad(idB):
        for i in range(srcB.propCount()):
            prop = srcB.propName(i)
            col = dst.propIndex(prop+'B')
            dst.setProp(id,col,srcB.getProp(idB,i))
    return id

def make_alchemical(atable, btable, ctable, block, keeper=None):
    ''' merge forcefield entries from m2 into m1. '''
    if btable is None: return

    b_constrained = btable.termPropIndex('constrained')
    c_constrained = None
    if not _msys.bad(b_constrained):
        c_constrained = ctable.addTermProp('constrained', int)
    params=ctable.params()
    tableBname = "alchemical_" + ctable.name()
    tableB = ctable.system().table(tableBname)
    if tableB is None:
        tableB = ctable.system().addTableFromSchema(tableBname,tableBname)
    paramsB = tableB.params()
                                    
    for (ta,tb), atoms in block:
        if ta>0 and tb==-1:
            # Keep the A state intact, no alchemical transformation
            pass

        elif ta>0 and tb==0:
            '''disappear the B state '''
            p = copy_alchemical(paramsB,
                                params,  ctable.param(ta-1),
                                None,    _msys.BadId)
            if keeper is not None:
                colA = paramsB.propIndex(keeper + 'A')
                colB = paramsB.propIndex(keeper + 'B')
                paramsB.setProp(p,colB,paramsB.getProp(p,colA))
            tableB.addTerm(ctable.atoms(ta-1), p)
            ctable.delTerm(ta-1)

        elif ta==-1 and tb>0:
            ''' copy parameters from state B '''
            ids=_msys.IdList()
            for a in atoms: ids.append(a)
            param=copy_param(params, btable.params(), btable.param(tb-1))
            term = ctable.addTerm( ids, param )
            if c_constrained is not None:
                cons = btable.getTermProp(tb-1, b_constrained)
                ctable.setTermProp(term, c_constrained, cons)

        elif ta==0 and tb>0:
            ''' disappear the A state '''
            p = copy_alchemical(paramsB,
                                None, _msys.BadId,
                                btable.params(), btable.param(tb-1))

            if keeper is not None:
                colA = paramsB.propIndex(keeper + 'A')
                colB = paramsB.propIndex(keeper + 'B')
                paramsB.setProp(p,colA,paramsB.getProp(p,colB))

            ids=_msys.IdList()
            for a in atoms: ids.append(a)
            tableB.addTerm(ids, p)

        elif ta>0 and tb>0:
            ''' A state morphs to B state '''
            # create a new entry in paramsB with values taken from A and B
            #paramB = copy_param(params, btable.params(), btable.param(tb-1))
            p = copy_alchemical(paramsB, 
                                params,          ctable.param(ta-1),
                                btable.params(), btable.param(tb-1))
            #ctable.setParamB(term, paramB)
            tableB.addTerm(ctable.atoms(ta-1), p)
            ctable.delTerm(ta-1)

        else:
            raise ValueError, "Unsupported mapping in %s: ta=%d, tb=%d" % (
                    ctable.name(), ta, tb)

        #if has_constrained and item.entryA and item.entryB:
            #item["constrained"] = 0


def make_alchemical_excl(table, block):
    for (ta,tb), atoms in block:
        ids=_msys.IdList()
        for i in sorted(atoms): ids.append(i)
        if ta>0 and (tb>0 or tb==-1):
            pass # already have this exclusion
        elif ta==-1 and (tb>0 or tb==-1):
            table.addTerm(ids, _msys.BadId)
    
def make_constraint_map( bmap, C, B):
    ''' provide just the alchemical part of m1 and m2.  Yields a lazy 
    sequence of corresponding items from m1 and m2; the m1 item may be
    None. '''

    # hash all the constraints in C from all categories by first atom
    consmap=dict()
    for name in C.tableNames():
        atable = C.table(name)
        if atable.category!=_msys.Category.constraint: continue
        for t in atable.terms():
            a=atable.atom(t,0)
            if a in consmap:
                raise ValueError, "Overlapping constraint for atom %d" % a
            consmap[a]=(name, t)

    # find constraints in m2 that overlap with m1
    for name in B.tableNames():
        btable = B.table(name)
        if btable.category!=_msys.Category.constraint: continue
        atable=C.table(name)
        for t in btable.terms():
            yield consmap.get( bmap[btable.atom(t,0)]), (name, t)

def create_ahn_constraint(mol, n):
    name='constraint_ah%d' % n
    return mol.addTableFromSchema(name, name)

def identical_params( params, i, j):
    ''' do the params in rows i and j have the same values? '''
    if i==j: return True
    if _msys.bad(i) or _msys.bad(j):
        return False
    for n in range(params.propCount()):
        if params.getProp(i,n) != params.getProp(j,n):
            return False
    return True

def MakeAlchemical(A, B, pairs):

    # put all the dummy atoms for A at the end
    # TODO: having made this transformation on pairs, we can probably 
    # simplify a fair bit of other things in this routine.
    pairs.sort()
    zero=None
    for i, (ai,aj) in enumerate(pairs):
        if ai==0:
            break
    pairs.extend(pairs[:i])
    pairs=pairs[i:]

    nC = len(pairs)
    apairs, bpairs = zip(*pairs)

    # clone just the alchemical parts of A and B
    atoms = _msys.IdList()
    for a in sorted(apairs):
        if a >= 0:
            atoms.append(a)
    C=_msys.Clone(A, atoms)

    batoms = _msys.IdList()
    for b in sorted(bpairs):
        if b >= 0:
            batoms.append(b)
    B=_msys.Clone(B, batoms)


    # add custom atom properties from B
    for i in range(B.atomPropCount()):
        C.addAtomProp(B.atomPropName(i), B.atomPropType(i))
    
    amap = _msys.IdList()
    bmap=dict()

    nbB = B.table('nonbonded')
    nbC = C.table('nonbonded')
    Czero = None
    alc = C.addTable('alchemical_nonbonded', 1, nbC.params())
    alc.category = _msys.parse_category('nonbonded')
    alc.addTermProp('chargeB', float)
    alc.addTermProp('moiety', int)

    for i, (ai, bi) in enumerate(pairs):
        assert ai>=0 or bi>=0
        if ai<0:
            # find the chain and residue in the corresponding B atom
            batm = B.atom(bi)
            bres = B.residue(batm.residue)
            bchn = B.chain(bres.chain)
            chn, res = find_residue(C, bres.name, bres.resid, bchn.name)
            if chn<0: 
                chn = C.addChain()
                C.chain(chn).name = bchn.name
            if res<0: 
                res = C.addResidue(chn)
                C.residue(res).name = bres.name
                C.residue(res).resid = bres.resid
            atm = C.addAtom(res)
            atom = C.atom(atm)
            atom.name = batm.name
            atom.mass = batm.mass
            atom.x = batm.x
            atom.y = batm.y
            atom.z = batm.z
            atom.vx = batm.vx
            atom.vy = batm.vy
            atom.vz = batm.vz
            atom.atomic_number = batm.atomic_number
            for p in range(B.atomPropCount()):
                prop = B.atomPropName(p)
                C.setAtomProp(atm, prop, B.getAtomProp(bi, prop))
            amap.append(atm)

            # make a zero term for the nonbonded
            if Czero is None: Czero = nbC.params().addParam()
            t = nbC.addTerm([atm], Czero)
            p = copy_param(nbC.params(), nbB.params(), nbB.param(bi))
            alcid = alc.addTerm([atm], p)
            alc.setTermProp(alcid,0,batm.charge)
        else:
            atom = C.atom(ai)
            amap.append(ai)
            if bi>=0:
                batm = B.atom(bi)
                p = copy_param(nbC.params(), nbB.params(), nbB.param(bi))
                alcid = alc.addTerm([ai], p)
                alc.setTermProp(alcid,0,batm.charge)
            else:
                if Czero is None: Czero = nbC.params().addParam()
                alcid = alc.addTerm([ai], Czero)
        if bi>=0:
            bmap[bi] = i

    # at this point we have all the dummies.  Now clone C again to put
    # them in the right order.
    C = _msys.Clone(C, amap)

    # stretches, angles, dihedrals
    ff='stretch_harm'
    bondmaps = make_block(C,amap, bmap, apairs, bpairs, C.table(ff), B.table(ff))
    ff='angle_harm'
    anglmaps = make_block(C,amap, bmap, apairs, bpairs, C.table(ff), B.table(ff))
    ff='dihedral_trig'
    dihemaps = make_block(C,amap, bmap, apairs, bpairs, C.table(ff), B.table(ff))
    ff='pair_12_6_es'
    pairmaps = make_pairmaps(C,bmap, C.table(ff), B.table(ff))

    kept.stage2( C, B, pairs, bmap, bondmaps, anglmaps, dihemaps, False )

    ff='stretch_harm'
    make_alchemical(A.table(ff), B.table(ff), C.table(ff), bondmaps, "r0")
    ff='angle_harm'
    make_alchemical(A.table(ff), B.table(ff), C.table(ff), anglmaps)
    ff='dihedral_trig'
    make_alchemical(A.table(ff), B.table(ff), C.table(ff), dihemaps)
    ff='pair_12_6_es'
    make_alchemical(A.table(ff), B.table(ff), C.table(ff), pairmaps)

    # bonds
    for b in B.bonds():
        bnd = B.bond(b)
        bi, bj = bmap[bnd.i], bmap[bnd.j]
        C.addBond(bi,bj)

    # exclusions and pairs
    ff='exclusion'
    exclmaps = make_exclmap(bmap, apairs, bpairs, C.table(ff), B.table(ff))
    exclmaps.extend(add_more_exclusions( pairs, bmap, A, B ))
    make_alchemical_excl(C.table(ff), exclmaps)

    # constraints
    for ca, (bname, tb) in make_constraint_map(bmap, C, B):
        btable = B.table(bname)
        bparams = btable.params()
        bp = btable.param(tb)
        ids=_msys.IdList()
        for a in btable.atoms(tb): ids.append(bmap[a])
        if ca is None:
            # copy constraint from B to C
            atable = create_ahn_constraint(C, len(ids)-1)
            param = copy_param(atable.params(), bparams, btable.param(tb))
            atable.addTerm( ids, param )

        else:
            # merge constraint from B into C
            cname, tc = ca
            ctable = C.table(cname)
            cparams = ctable.params()
            nparams = cparams.propCount()
            p = ctable.param(tc)
            cids = ctable.atoms(tc)
            props = [cparams.getProp(p,i) for i in range(nparams)]
            for i,b in enumerate(ids[1:]):
                if b not in cids:
                    cids.append(b)
                    props.append(bparams.getProp(bp, i))
            ctable.delTerm(tc)
            mtable = create_ahn_constraint(C, len(props))
            mparams = mtable.params()
            mp = mparams.addParam()
            for i in range(len(props)):
                mparams.setProp(mp, i, props[i])
            mtable.addTerm(cids, mp)


    # tack on the non-alchemical part
    alc_ids = A.atoms()[len(atoms):]
    alc = _msys.Clone(A, alc_ids)
    C.append(alc)

    # coalesce and clone to obtain the minimal set of params
    C.coalesceTables()
    return _msys.Clone(C, C.atoms())

