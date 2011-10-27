
'''
fep: generate a new ent with alchemical terms.

'''

import kept
from math import sqrt
import sys

'''
support for generation of full pairs
'''

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

def find_nonbonded_entry(nb, atom):
    for item in nb.items:
        if item.atom(0)==atom:
            return item.entry
    raise ValueError, "Missing nonbonded term for atom %s" % atom

def configure_vdw_12_6_entry(entry, vdw_rule, entry_i, entry_j):
    ''' Using the given nonbonded entries, configure the vdw part of
    the given pair entry. '''
    # get sigma, epsilon for each atom
    props=('sigma', 'epsilon')
    param_i = [entry_i[p] for p in props]
    param_j = [entry_j[p] for p in props]
    if vdw_rule=='geometric':
        aij, bij = combine_geometric(param_i, param_j)
    elif vdw_rule=='arithmetic/geometric':
        aij, bij = combine_arith_geom(param_i, param_j)
    else:
        raise ValueError, "unsupported rule '%s' for vdw_12_6" % vdw_rule
    entry['aij'] = aij
    entry['bij'] = bij

def configure_vdw_exp_6_entry(entry, vdw_rule, entry_i, entry_j):
        props=('alpha', 'epsilon', 'rmin')
        param_i = [entry_i[p] for p in props]
        param_j = [entry_j[p] for p in props]
        vdw_rule = nb['vdw_rule'].lower()
        if vdw_rule=='lb/geometric':
            aij, bij, cij = combine_lb_geometric(param_i, param_j)
        else:
            raise ValueError, "unsupported rule '%s' for vdw_exp_6" % vdw_rule
        entry['aij'] = aij
        entry['bij'] = bij
        entry['cij'] = cij

def get_pairs_table(mol):
    nb = mol.getNBody('nonbonded')
    vdw_funct = nb['vdw_funct'].lower()
    if vdw_funct=='vdw_12_6':
        pname = 'pair_12_6_es'
    elif vdw_funct=='vdw_exp_6':
        pname = 'pair_exp_6_es'
    else:
        raise ValueError, "Unsupported vdw_funct '%s'" % vdw_funct

    pairs = mol.getNBody(pname)
    if not pairs:
        raise ValueError, "No pair_12_6_es table, have %s" % ( 
                mol.getNBodyNames())
    return pairs

def make_full_pair_entry(mol, iatom, jatom):
    ''' if an item already exists for these atoms, use that item, otherwise
    create a new one with the full entry in the entryA slot.
    '''
    pairs = get_pairs_table(mol)
    entry = pairs.addEntry()
    entry['qij'] = iatom.charge * jatom.charge

    nb = mol.getNBody('nonbonded')
    nb_i=find_nonbonded_entry(nb, iatom)
    nb_j=find_nonbonded_entry(nb, jatom)
    vdw_funct = nb['vdw_funct'].lower()
    vdw_rule = nb['vdw_rule'].lower()

    if vdw_funct=='vdw_12_6':
        configure_vdw_12_6_entry( entry, vdw_rule, nb_i, nb_j )
    elif vdw_funct=='vdw_exp_6':
        configure_vdw_exp_6_entry( entry, vdw_rule, nb_i, nb_j)
    else:
        raise_ValueError, "Unsupported vdw_funct %s'" % vdw_funct

    return pairs, entry

def canonical(ids):
    n=len(ids)
    if   n==2: a,b=ids
    elif n==3: a,b=ids[0],ids[2]
    elif n==4: a,b=ids[1],ids[2]
    else: raise ValueError, "Expect len(ids) in (2,3,4), got %d" % len(ids)
    if a>b:
        ids.reverse()
    return ids

def keep_dummies(imap, bmap):
    '''Assign kept (-1) and not kept (0) for terms involving only dummy atoms.
    '''
    for inds, term in bmap:
        ti, tj = inds
        if tj<1:
            dums=[imap[x] for x in term if imap[x]<0]
            inds[1] = -1 if len(dums)==len(term) else 0
        elif ti<1:
            dums=[x for x in term if x<0]
            inds[0] = -1 if len(dums)==len(term) else 0

def make_block( map, imap, nbA, nbB ):
    # a_item_list: the ids of the atoms in each atom, canonicalized
    natoms = nbA.atoms_per_item
    arange = range(natoms)
    block=[]

    # bterms: mapping from canonicalized atoms to 1-based index
    bterms={}
    for i, item in enumerate(nbB.items):
        ids=[item.atom(j).id+1 for j in arange]
        key=canonical([map[x] for x in ids])
        bterms[tuple(key)]=i+1

    for i, item in enumerate(nbA.items):
        key=canonical([item.atom(j).id+1 for j in arange])
        j = bterms.pop(tuple(key),0)
        block.append( [[i+1, j], key] )

    # add remaining items in b_item_dict
    block.extend( sorted( [[0,j], list(ids)] for ids, j in bterms.items() ))

    keep_dummies( imap, block )

    return block

def make_pairmaps(map, nbA, nbB):
    # make a and b items same way as exclmap
    a_item_list=[sorted(
        (item.atom(0).id+1, item.atom(1).id+1)) for item in nbA.items]
    b_item_dict=dict()
    for i,x in enumerate(nbB.items):
        item=sorted( (map[x.atom(0).id+1], map[x.atom(1).id+1]) )
        b_item_dict[tuple(item)]=i+1

    block=list()
    for i, item in enumerate(a_item_list):
        j=b_item_dict.get(tuple(item), 0) # 0 instead of -1 as in exclusion
        block.append( [[i+1,j], item] )
        if j!=0:
            del b_item_dict[tuple(item)]

    block2=list()
    for item, j in b_item_dict.items():
        block2.append( [[0,j], item] )
    block2.sort()
    block.extend(block2)

    return block

def make_exclmap( map, imap, invmap, exclA, exclB, entA, entB ):
    a_item_list=[sorted(
        (item.atom(0).id+1, item.atom(1).id+1)) for item in exclA.items]
    b_item_dict=dict()
    for i,x in enumerate(exclB.items):
        item=sorted( (map[x.atom(0).id+1], map[x.atom(1).id+1]) )
        b_item_dict[tuple(item)]=i+1

    Aatoms=[a for a in entA.atoms]
    Batoms=[a for a in entB.atoms]

    block=list()
    for i, item in enumerate(a_item_list):
        j=b_item_dict.get(tuple(item), -1)
        block.append( [[i+1,j], item] )
        if j!=-1:
            del b_item_dict[tuple(item)]
        else:
            # An excluded pair of atoms in A mapping onto an unexcluded pair
            # of atoms in B: we keep the exclusion, but add an artificial
            # pair to B.
            ai=imap[item[0]]
            aj=imap[item[1]]
            if item[0]>0 and item[1]>0 and ai>0 and aj>0:
                #print 'Offset exclusion %d-%d in A by pair interaction in B: %d-%d' % (item[0], item[1], ai, aj)
                i_atom=Batoms[ai-1]
                j_atom=Batoms[aj-1]
                pairs, entry = make_full_pair_entry(entB, i_atom, j_atom)
                # we want the pairs table for the A state
                pairs = get_pairs_table(entA)
                for t in pairs.items:
                    if t.atom(0).id==item[0]-1 and t.atom(1).id==item[1]-1:
                        t.entryB = copy_entry(pairs, entry)
                        break
                else:
                    pairs.addItem( 
                            pairs.addEntry(), 
                            copy_entry(pairs, entry),
                            Aatoms[item[0]-1],
                            Aatoms[item[1]-1])

    # Add the remaining exclusions in B
    block2=list()
    for item, j in b_item_dict.items():
        block2.append( [[-1,j], item] )
        ai=imap[item[0]]
        aj=imap[item[1]]
        if item[0]>0 and item[1]>0 and ai>0 and aj>0:
            #print 'Offset exclusion %d-%d in B by pair interaction in A: %d-%d' % (ai, aj, item[0], item[1])
            i_atom=Aatoms[item[0]-1]
            j_atom=Aatoms[item[1]-1]
            pairs, entry = make_full_pair_entry(entA, i_atom, j_atom)
            for item in pairs.items:
                if item.atom(0)==i_atom and item.atom(1)==j_atom:
                    item.entry = entry
                    break
            else:
                pairs.addItem( entry, pairs.addEntry(), i_atom, j_atom )

    block2.sort()
    block.extend(block2)

    return block

def add_more_exclusions( block, pairs, imap, m1, m2, cutoff=6.0 ):
    ai_list=list()
    aj_list=list()
    atoms_1=[x for x in m1.atoms]
    atoms_2=[x for x in m2.atoms]
    cut2=cutoff*cutoff
    for ai,aj in pairs:
        if ai>0 and aj<0: ai_list.append(ai)
        if ai<0 and aj>0: aj_list.append(ai) # yes, ai; using imap later
    for ai in ai_list:
        ipos=atoms_1[ai-1].pos
        for aj in aj_list:
            jpos=atoms_2[imap[aj]-1].pos
            jpos-=ipos
            if geom.Length2(jpos) < cut2:
                block.append( [[-1,-1], [ai,aj]] )

def get_nbody(e, name):
    if name in e.getNBodyNames():
        return e.getNBody(name)
    return ent.NBody()

def alchemical_combine( m1, m2, pairs ):
    ''' merge m2 into m1.  Return mapping from aj into ai '''
    invmap = dict()
    atom_1=[x for x in m1.atomList()]
    atom_2=[x for x in m2.atomList()]
    nb1=get_nbody(m1,'nonbonded')
    nb2=get_nbody(m2,'nonbonded')
    # make a mapping from atoms to their nonbonded item
    nbitem1=dict()
    nbitem2=dict()
    for x in nb1.items: nbitem1[x.atom(0)] = x
    for x in nb2.items: nbitem2[x.atom(0)] = x

    ailist=[x for x,y in pairs]
    # check that we have a valid atom list
    if sorted(list(set(map(abs,ailist))))!=range(1,1+len(ailist)):
        raise ValueError, "Bad atom map, first column has duplicates"

    # copy atom props from m2 to m1
    props1 = set(m1.getAtomProps())
    props2 = set(m2.getAtomProps())
    newprops = props2.difference(props1)

    # FIXME: awaiting a better ent interface for this.
    for p in newprops:
        ptype = type(atom_2[0][p])
        if ptype is int: 
            dtype=ent.DICT_VALUE_INT
        elif ptype is float: 
            dtype=ent.DICT_VALUE_FLOAT
        else:
            dtype=ent.DICT_VALUE_STRING
        m1.addAtomProp(p, dtype)


    # make room for the dummy atoms
    ndummy = len(filter(lambda x: x<0, ailist))
    nreal = len(ailist)-ndummy
    for i in xrange(nreal,len(atom_1)):
        atom_1[i].id += ndummy

    for ai, aj in pairs:
        if ai>0 and aj>0:
            a1=atom_1[ai-1]
            a2=atom_2[aj-1]
            invmap[aj-1]=ai-1
            a1.alch=True
            a1.chargeB=a2.charge
            nbitem1[a1].entryB=copy_entry(nb1, nbitem2[a2])

        elif ai>0 and aj<0:
            a1=atom_1[ai-1]
            a1.alch=True
            nbitem1[a1].entryB=nb1.addEntry()

        elif ai<0 and aj>0:
            a2=atom_2[aj-1]
            invmap[aj-1]=abs(ai)-1
            # add a new dummy atom.  If possible, use the chain and residue 
            # from aj; otherwise create new ones.
            chain=m1.findChain(a2.chain.name)
            if not chain:
                chain=m1.addChain(a2.chain.name)
            residue=chain.findResidue(a2.residue.name, a2.residue.num)
            if not residue:
                residue=chain.appendResidue()
                residue.name=a2.residue.name
                residue.num=a2.residue.num
            a1=residue.appendAtom()
            a1.id=abs(ai)-1
            a1.alch=True
            a1.charge=0
            a1.chargeB=a2.charge
            a1.mass=a2.mass
            a1.pos=a2.pos
            a1.vel=a2.vel
            a1.name=a2.name
            a1.num=a2.num
            for p in newprops:
                a1[p]=a2[p]
            nb1.addItem(
                    nb1.addEntry(),
                    copy_entry( nb1, nbitem2[a2] ),
                    a1 )
        else:
            raiseValueError, "Illegal pair: ai=%d, aj=%d" % (ai,aj)

    m1.sortAtomsByID()
    #for k,v in invmap.items(): print "INV %d -> %d" % (k+1,v+1)
    return invmap

def fixup_bonds(m1, m2, invmap):
    ''' map bonds in m2 into m1 '''
    atom_1=[x for x in m1.atomList()]
    for a2 in m2.atoms:
        a_id=invmap.get(a2.id)
        if a_id is None: continue
        a1=atom_1[a_id]
        for b2 in a2.bondedAtoms():
            b_id=invmap[b2.id]
            b1=atom_1[b_id]
            if a1.id < b1.id:
                m1.addBond(a1,b1)
            else:
                m1.addBond(b1,a1)

def copy_entry( nb, entry ):
    ''' create a new entry in nbody, and copy properties from entry to
    the new entry for each of the properties in nbody. '''
    dst = nb.addEntry()
    for p in nb.prop_names: 
        val=entry[p]
        if val is None: continue
        dst[p]=val
    return dst

def make_alchemical( m1, m2, nbname, block, keeper=None ):
    ''' merge forcefield entries from m2 into m1. '''
    nb1=get_nbody(m1, nbname)
    nb2=get_nbody(m2, nbname)
    nb1.addItemProp('moiety', 0)
    atoms_1=[x for x in m1.atoms]
    items_1=[x for x in nb1.items]
    items_2=[x for x in nb2.items]
    has_constrained = nb2.hasItemProp("constrained")
    for (ta,tb), atoms in block:
        if ta>0 and tb==-1:
            # Keep the A state intact, no alchemical transformation
            pass

        elif ta>0 and tb==0:
            '''disappear the B state '''
            item=items_1[ta-1]
            item.entryB=copy_entry(nb1, item.entryA)
            for p in nb1.prop_names:
                if p!=keeper:
                    item.entryB[p]=0.0

        elif ta==-1 and tb>0:
            # copy parameters from state B
            item=nb1.addItem( 
                    copy_entry( nb1, items_2[tb-1].entry ),
                    *[atoms_1[abs(x)-1] for x in atoms] )
            if has_constrained:
                cons = items_2[tb-1]["constrained"]
                item["constrained"] = cons

        elif ta==0 and tb>0:
            ''' disappear the A state '''
            entryB=copy_entry( nb1, items_2[tb-1].entry )
            entryA=copy_entry( nb1, entryB )
            for p in nb1.prop_names:
                if p!=keeper:
                    entryA[p]=0.0
            item = nb1.addItem( 
                    entryA, entryB, *[atoms_1[abs(x)-1] for x in atoms] )

        elif ta>0 and tb>0:
            ''' A state morphs to B state '''
            item = items_1[ta-1]
            item.entryB = copy_entry( nb1, items_2[tb-1].entry )

        else:
            raise ValueError, "Unsupported mapping in %s: ta=%d, tb=%d" % (
                    nbname, ta, tb)

        if has_constrained and item.entryA and item.entryB:
            item["constrained"] = 0


def make_alchemical_exclusions(m,nb,block):
    # need a dummy entry for the exclusion items
    entry=nb.addEntry()
    atoms_1=[a for a in m.atoms]
    for (ta,tb), atoms in block:
        ai,aj=sorted([abs(a)-1 for a in atoms])
        if ta>0 and (tb>0 or tb==-1):
            pass # already have this exclusion
        elif ta==-1 and (tb>0 or tb==-1):
            nb.addItem(entry, atoms_1[ai], atoms_1[aj])
        else:
            raise ValueError, "Unsupported exclusion: ta=%d, tb=%d" % (ta,tb)


def make_constraint_map( invmap, m1, m2 ):
    ''' provide just the alchemical part of m1 and m2.  Yields a lazy 
    sequence of corresponding items from m1 and m2; the m1 item may be
    None. '''

    # hash all the constraints in m1 from all categories by first atom
    consmap=dict()
    for nbname in m1.getNBodyNames():
        nb1=m1.getNBody(nbname)
        if nb1.category!='constraint': continue
        for item in nb1.items:
            aid=item.atom(0).id
            if aid in consmap:
                raise ValueError, "Found overlapping constraint for atom %d" % aid
            consmap[aid]=(nbname, item)

    # find constraints in m2 that overlap with m1
    for nbname in m2.getNBodyNames():
        nb2=m2.getNBody(nbname)
        natoms=nb2.atoms_per_item
        if nb2.category!='constraint': continue
        nb1=m2.getNBody(nbname)
        assert nb1, "FIXME: Constraint type %s in B state but not in A state" % nbname
        for itemB in nb2.items:
            yield consmap.get( invmap[itemB.atom(0).id] ), (nbname, itemB)

def create_ahn_constraint(mol, n):
    name='constraint_ah%d' % n
    nb=mol.getNBody(name)
    if not nb:
        nb=mol.createNBody(name, n+1)
        nb.category='constraint'
        for i in range(1,n+1):
            nb.addProp('r%d' % i, 0.0)
        nb.addItemProp('moiety', 0)
    return nb

def MakeAlchemical(entA, entB, pairs, keep_extra_dihedral=False):

    map, imap = kept.create_map(pairs)
    
    # Construct views for the alchemical parts of the two entities
    selA='index ' + ' '.join([str(x-1) for x in imap.keys() if x>0])
    selB='index ' + ' '.join([str(x-1) for x in  map.keys() if x>0])
    alcA=entA.select(selA)
    alcB=entB.select(selB)

    ff = 'stretch_harm'
    block_bondmaps = make_block(map,imap,get_nbody(alcA,ff),get_nbody(alcB,ff))

    ff = 'angle_harm'
    block_anglemaps = make_block(map,imap,get_nbody(alcA,ff),get_nbody(alcB,ff))

    ff = 'dihedral_trig'
    block_dihedmaps = make_block(map,imap,get_nbody(alcA,ff),get_nbody(alcB,ff))

    # update the blocks using the kept terms algorithm
    kept.stage2( alcA, alcB, pairs,
                 block_bondmaps, block_anglemaps, block_dihedmaps,
                 keep_extra_dihedral)

    # at this point the bond, angle, and dihedral maps are finished.

    # particles
    invmap = alchemical_combine( entA, entB, pairs )

    # bonds
    fixup_bonds(entA, entB, invmap)

    # stretch
    make_alchemical(entA, entB, 'stretch_harm', block_bondmaps, "r0")

    # angle
    make_alchemical(entA, entB, 'angle_harm', block_anglemaps )

    # dihedral
    make_alchemical(entA, entB, 'dihedral_trig', block_dihedmaps)

    # pairs
    ff = 'pair_12_6_es'
    block_pairmaps = make_pairmaps(
            map, 
            get_nbody(alcA,ff), 
            get_nbody(alcB,ff))
    make_alchemical(entA, entB, ff, block_pairmaps)

    # exclusions and pairs
    ff = 'exclusion'
    block_excl = make_exclmap(
            map, imap, invmap,
            get_nbody(alcA,ff), get_nbody(alcB,ff), entA, entB)
    add_more_exclusions( block_excl, pairs, imap, alcA, alcB)
    make_alchemical_exclusions(entA, get_nbody(entA,ff), block_excl)

    # constraints
    Aatoms=[a for a in entA.atoms]
    for A, (nbnameB, itemB) in make_constraint_map(invmap, alcA, alcB):
        nbB=entB.getNBody(nbnameB)
        ids=[invmap[itemB.atom(i).id] for i in range(nbB.atoms_per_item)]
        if A is None:
            # copy itemB into nbA
            nbA=create_ahn_constraint(entA,len(ids)-1)
            nbA.addItem(
                copy_entry(nbA, itemB),
                *[Aatoms[i] for i in ids])
        else:
            # merge itemB into itemA
            nbnameA, itemA = A
            nbA=entA.getNBody(nbnameA)
            aids=[itemA.atom(i).id for i in range(nbA.atoms_per_item)]
            props=[itemA['r%d' % i] for i in range(1,nbA.atoms_per_item)]
            for i, b in enumerate(ids[1:]):
                if b not in aids:
                    aids.append(b)
                    props.append(itemB['r%d' % (i+1)])
            #print "Creating merged constraint", aids, props
            nbA.deleteItem(itemA)
            nbB.deleteItem(itemB)
            nb=create_ahn_constraint(entA, len(props))
            entry=nb.addEntry()
            for i in range(len(props)):
                entry['r%d' % (i+1)] = props[i]
            nb.addItem(entry, *[Aatoms[i] for i in aids])

def UnchargeAlchemical( mol, seltext ):
    ''' make the selected atoms alchemical, and set their charge to 0 in the
    B state.
    '''
    nonbonds=mol.table('nonbonded')
    nbmap=dict()
    if nonbonds:
        for term in nonbonds.terms:
            nbmap[term.atom(0)] = term

    atoms=mol.atomselect(seltext)
    if not atoms:
        print >> sys.stderr, "WARNING: UnchargeAlchemical -- no atoms in selection"
        return
    for atom in atoms:
        atom.alchemical=True
        atom.chargeB = 0.0
        nb = nbmap.get(atom)
        if nb is not None:
            nb.paramB = nb.param
        elif nonbonds is not None:
            raise ValueError, "No nonbond found for atom %s" % atom


