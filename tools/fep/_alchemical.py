
import _msys
from msys import kept

def find_residue(mol, resname, resid, chainname):
    for chn in mol.chains():
        if mol.chain(chn).name == chainname:
            for res in mol.residuesForChain(chn):
                residue = mol.residue(res)
                if residue.name == resname and residue.num == resid:
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

def make_block( amap, bmap, apairs, bpairs, atable, btable):
    # a_item_list: the ids of the atoms in each atom, canonicalized
    natoms = atable.atomCount()
    arange = range(natoms)
    block=[]

    # bterms: mapping from canonicalized atoms to 1-based index
    bterms={}
    for i in range(btable.termCount()):
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

def make_pairmaps(bmap, atable, btable):
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

def make_alchemical(atable, btable, ctable, block, keeper=None):
    ''' merge forcefield entries from m2 into m1. '''
    b_constrained = btable.termPropIndex('constrained')
    c_constrained = None
    if not _msys.bad(b_constrained):
        c_constrained = ctable.addTermProp('constrained', int)
    params=ctable.params()
    for (ta,tb), atoms in block:
        if ta>0 and tb==-1:
            # Keep the A state intact, no alchemical transformation
            pass

        elif ta>0 and tb==0:
            '''disappear the B state '''
            term=ta-1
            paramB=copy_param(params, params, ctable.param(term))
            ctable.setParamB(term, paramB)
            for index in range(params.propCount()):
                if params.propName(index)!=keeper:
                    params.setProp(paramB, index, 0.0)

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
            ids=_msys.IdList()
            for a in atoms: ids.append(a)
            param=copy_param(params, btable.params(), btable.param(tb-1))
            paramB=params.duplicate(param)
            for index in range(params.propCount()):
                if params.propName(index)!=keeper:
                    params.setProp(param, index, 0.0)
            term = ctable.addTerm( ids, param )
            ctable.setParamB(term, paramB)

        elif ta>0 and tb>0:
            ''' A state morphs to B state '''
            term=ta-1
            paramB = copy_param(params, btable.params(), btable.param(tb-1))
            ctable.setParamB(term, paramB)

        else:
            raise ValueError, "Unsupported mapping in %s: ta=%d, tb=%d" % (
                    ctable.name(), ta, tb)

        #if has_constrained and item.entryA and item.entryB:
            #item["constrained"] = 0

def MakeAlchemical(A, B, pairs):

    nC = len(pairs)
    apairs, bpairs = zip(*pairs)

    # clone just the alchemical part of A
    atoms = _msys.IdList()
    for a in apairs:
        if a >= 0:
            atoms.append(a)
    C=_msys.Clone(A, atoms)
    amap = _msys.IdList()
    bmap=dict()

    nbB = B.table('nonbonded') if 'nonbonded' in B.tableNames() else None
    nbC = C.table('nonbonded') if 'nonbonded' in C.tableNames() else None
    Czero = None

    for i, (ai, bi) in enumerate(pairs):
        assert ai>=0 or bi>=0
        if ai<0:
            # find the chain and residue in the corresponding B atom
            batm = B.atom(bi)
            bres = B.residue(batm.residue)
            bchn = B.chain(bres.chain)
            chn, res = find_residue(C, bres.name, bres.num, bchn.name)
            if chn<0: chn = C.addChain()
            if res<0: res = C.addResidue(chn)
            atm = C.addAtom(res)
            atom = C.atom(atm)
            atom.alchemical = True
            atom.mass = batm.mass
            atom.x = batm.x
            atom.y = batm.y
            atom.z = batm.z
            atom.chargeB = batm.charge
            atom.atomic_number = batm.atomic_number
            amap.append(atm)
            if nbC is not None:
                assert nbB is not None
                # make a zero term for the nonbonded
                if Czero is None: Czero = nbC.params().addParam()
                t = nbC.addTerm([atm], Czero)
                p = copy_param(nbC.params(), nbB.params(), nbB.param(bi))
                nbC.setParamB(t, p)
        else:
            atom = C.atom(ai)
            amap.append(ai)
            atom.alchemical = True
            if bi>=0:
                batm = B.atom(bi)
                atom.chargeB = batm.charge
                p = copy_param(nbC.params(), nbB.params(), nbB.param(bi))
                nbC.setParamB(ai, p)
            else:
                if Czero is None: Czero = nbC.params().addParam()
                nbC.setParamB(ai, Czero)
        if bi>=0:
            bmap[bi] = i

    # at this point we have all the dummies.  Now clone C again to put
    # them in the right order.
    C = _msys.Clone(C, amap)

    # stretches, angles, dihedrals
    ff='stretch_harm'
    bondmaps = make_block(amap, bmap, apairs, bpairs, C.table(ff), B.table(ff))
    ff='angle_harm'
    anglmaps = make_block(amap, bmap, apairs, bpairs, C.table(ff), B.table(ff))
    ff='dihedral_trig'
    dihemaps = make_block(amap, bmap, apairs, bpairs, C.table(ff), B.table(ff))
    ff='pair_12_6_es'
    pairmaps = make_pairmaps(bmap, C.table(ff), B.table(ff))

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
        if apairs[bi]<0 or apairs[bj]<0:
            C.addBond(bi,bj)

    # tack on the non-alchemical part
    alc_ids = A.atoms()[len(atoms):]
    alc = _msys.Clone(A, alc_ids)
    C.append(alc)

    # reassign gids in the order they're given by the input maps 
    for a in C.atoms(): C.atom(a).gid = a

    return C

