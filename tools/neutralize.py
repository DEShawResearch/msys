from __future__ import print_function

import msys, math
import random

def compute_center(residue):
    tm=0.0
    tx=0.0
    ty=0.0
    tz=0.0
    for a in residue.atoms:
        m = a.mass
        tm += m
        tx += m*a.x
        ty += m*a.y
        tz += m*a.z
    if tm:
        tx /= tm
        ty /= tm
        tz /= tm
    return (tx,ty,tz)

def dist2(pi, pj):
    d=0.0
    for i,j in zip(pi, pj):
        d += (i-j)**2
    return math.sqrt(d)

def parse_ion(name):
    if isinstance(name, str):
        anum = msys.ElementForAbbreviation(name.capitalize())
    elif isinstance(name, msys.System):
        return name
    else:
        anum = int(name)
    mol = msys.CreateSystem()
    atm = mol.addAtom()
    atm.atomic_number = anum
    atm.name = msys.AbbreviationForElement(anum)
    grp = msys.GroupForElement(anum)
    if grp == 1:
        atm.charge = 1
    elif grp == 2:
        atm.charge = 2
    elif grp == 17:
        atm.charge = -1
    else:
        raise RuntimeError("Could not determine for element '%s'" % name)
    atm.formal_charge = int(atm.charge)
    return mol

def Neutralize(mol, cation='Na', anion='Cl', 
        charge='formal_charge',
        chain='ION', chain2='ION2',
        solute_pad=5.0,
        ion_pad=3.0,
        concentration=0.0,
        keep='none',
        random_seed=0,
        verbose=False):

    """
    Replace water molecules with ions

    Arguments:
        mol (msys.System): input solvated system
        cation (str): element for cation species
        anion (str): element for anion species
        charge (str): either 'formal_charge', 'charge', or ZZ

        chain (str): chain name for counterions
        chain2 (str): chain name for counter-counterions
        concentration (float): 
        solute_pad (float): minimum distance between placed ions and non-water
        ion_pad (float): minimum distance between placed ions.
        concentration (float): molar concentration of added ions with same charge
        keep (str): selection of ions/waters that cannot be deleted or replaced
        random_seed (int): random seed to determine which water molecules to replace

    Returns:
        msys.System: a new system with ions added

    Notes:
        It is up to the caller to ensure that appropriate charges are
        assigned to the atoms in the solute, corresponding to the
        'charge' argument.  Use msys.AssignBondOrderAndFormalCharge if
        formal charges are missing.

        The number of ions with the same charge as the solute is calculated as

              nother = int((conc / 55.345) * (ntotalwat - nions))

        where ntotalwat is the total number of waters in the original
        system and nions is the number of counterions needed to achieve
        charge neutrality.
    """

    # Clone system
    mol = mol.clone()

    # create a single ct for all the ions
    ct = mol.addCt()
    ct.name = 'ion'

    # first, we add sufficient counterions to neutralize.  Then we add 
    # counterions and counter-counterions until counter-counterions are
    # up to the desired concentration.  
    cationsys = parse_ion(cation)
    anionsys = parse_ion(anion)
    cationatom = cationsys.atoms[0]
    anionatom = anionsys.atoms[0]

    solute=mol.select("""not ((atomicnumber %d and not bonded) or (atomicnumber %d and not bonded)) or (%s)""" %
                      (cationatom.atomic_number, anionatom.atomic_number, keep))
    if charge == 'formal_charge':
        cg = sum(a.formal_charge for a in solute)
    elif charge == 'charge':
        cg = sum(a.charge for a in solute)
    else:
        cg = float(charge)
    if verbose:
        print("neutralize: solute charge=%s" % cg)

    if cg >= 0:
        ionsys = anionsys
        othersys = cationsys
        ionatom = anionatom
        otheratom = cationatom
    else:
        ionsys = cationsys
        othersys = anionsys
        ionatom = cationatom
        otheratom = anionatom

    nions = int(math.fabs(cg/ionatom.charge)+0.5)

    # find the water residues
    water = mol.select('water and (not hydrogen) and (not within %f of (not water)) and (not (%s))'
            % (solute_pad, keep))
    residues = sorted(set(a.residue for a in water), key=lambda x: x.id)
    nwat = len(residues)
    if verbose:
        print("waters available to be replaced by ions:", nwat)

    # compute number of ions already present
    nions_prev = len(mol.select('atomicnumber %d and not bonded and (not (%s))'
                                % (ionatom.atomic_number, keep)))
    nother_prev = len(mol.select('atomicnumber %d and not bonded and (not (%s))'
                                 % (otheratom.atomic_number,keep)))

    if verbose:
        print("Starting with %d %s ions" % (nions_prev, ionatom.name))
        print("Starting with %d %s ions" % (nother_prev, otheratom.name))

    # convert molar concentration to number based on available waters.  The
    # molar concentration of water is about 55.345 mol/L.  Use all available
    # waters to calculate the number of ions to add.
    ntotalwat = len(set(a.residue for a in mol.select('water')))
    if verbose:
        print("Starting with %d water molecules" % ntotalwat)
    cgratio = math.fabs(otheratom.charge/ionatom.charge)
    nother = int((concentration/55.345) * (ntotalwat-nions+nions_prev))
    nions += int(nother*cgratio)


    if nions >= 0 and verbose:
        print("New system should contain %d %s ions" % (nions, ionatom.name))
    if nother >= 0 and verbose:
        print("New system should contain %d %s ions" % (nother, otheratom.name))

    # subtract off the ions already present in solution
    nions -= nions_prev
    nother -= nother_prev
    if nions < 0:
        # delete ions
        sel=mol.select('atomicnumber %d and not bonded and not (%s)' %
                (ionatom.atomic_number, keep))
        if len(sel) < -nions:
            raise RuntimeError("""Cannot decrease concentration - 
            not enough ions to delete""")
        for r in sel[:-nions]: r.remove()
        nions = 0
    if nother < 0:
        # delete other 
        sel=mol.select('atomicnumber %d and not bonded and not (%s)' %
                (otheratom.atomic_number, keep))
        if len(sel) < -nother:
            raise RuntimeError("""Cannot decrease concentration - not enough 
            other ions to delete""")
        for r in sel[:-nother]: r.remove()
        nother = 0


    if nwat < nions + nother:
        raise RuntimeError("""Only %d waters found; not enough to 
        neutralize""" % nwat)

    # Shuffle the residues
    random.seed(random_seed)
    random.shuffle(residues)

    # Remove overly close residues among the first nions + nother waters
    # Cache the centers to save time
    centers={}
    ionpad2 = ion_pad * ion_pad
    for i in range(nions + nother):
        ri = residues[i]
        try: pi = centers[ri.id]
        except KeyError:
            pi = centers[ri.id] = compute_center(ri)
        j = i+1
        while j < nions+nother:
            if len(residues) < nions+nother:
                raise RuntimeError("Not enough waters or too large ion_pad.")
            rj = residues[j]
            try: pj = centers[rj.id]
            except KeyError:
                pj = centers[rj] = compute_center(rj)
            d2 = dist2(pi,pj)
            if d2 < ionpad2:
                del residues[j]
            else:
                j += 1

    keep_atoms = set([x.id for x in mol.select('(%s)' % (keep,))])
    keep_residues = set([x.residue.id for x in mol.select('(%s)' % (keep,))])
    residues_removed = set()

    def multiple_waters(residue):
        real_atoms = [atom for atom in residue.atoms if atom.atomic_number > 0]
        return len(real_atoms) > 3

    def handle_multiple_waters(residue):
        raise RuntimeError("Some residues (%d) contain multiple waters. Use dms-fix-water-residues." % residue.id)

    if nions > 0:
        for i in range(nions):
            res=residues[i]
            if multiple_waters(res):
                handle_multiple_waters(res)
            newion = ct.append(ionsys)[0]
            newion.pos = compute_center(res)
            newion.residue.chain.name = chain
            newion.residue.resid=i+1
            residues_removed.add(res.id)
            res.remove()

    if nother > 0:
        for i in range(nions, nions+nother):
            res=residues[i]
            if multiple_waters(res):
                handle_multiple_waters(res)
            newion = ct.append(othersys)[0]
            newion.pos = compute_center(res)
            newion.residue.chain.name = chain2
            newion.residue.resid=i+1-nions
            residues_removed.add(res.id)
            res.remove()

    assert not (residues_removed & keep_residues)

    mol.coalesceTables()
    return mol.clone()

