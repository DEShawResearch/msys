
import _msys

def find_residue(mol, resname, resid, chainname):
    for chn in mol.chains():
        if mol.chain(chn).name == chainname:
            for res in mol.residuesForChain(chn):
                residue = mol.residue(res)
                if residue.name == resname and residue.num == resid:
                    return chn, res
            return chn, -1
    return -1, -1

def MakeAlchemical(A, B, pairs):

    nA = A.atomCount()
    nB = B.atomCount()
    nC = len(pairs)
    assert nA <= nC
    assert nB <= nC

    apairs, bpairs = zip(*pairs)

    C=_msys.Clone(A, A.atoms())
    amap = _msys.IdList()
    bmap=dict()
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
            atom.mass = batm.mass
            atom.x = batm.x
            atom.y = batm.y
            atom.z = batm.z
            atom.chargeB = batm.chargeB
            atom.atomic_number = batm.atomic_number
            amap.append(atm)
        else:
            amap.append(ai)
        if bi>=0:
            bmap[bi] = i

    # at this point we have all the dummies.  Now clone C again to put
    # them in the right order.
    C = _msys.Clone(C, amap)

    # bonds
    for b in B.bonds():
        bnd = B.bond(b)
        bi, bj = bmap[bnd.i], bmap[bnd.j]
        if apairs[bi]<0 or apairs[bj]<0:
            C.addBond(bi,bj)

    # reassign gids in the order they're given by the input maps 
    for a in C.atoms(): C.atom(a).gid = a

    return C

