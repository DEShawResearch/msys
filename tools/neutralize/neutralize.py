
import os, msys, math
from collections import namedtuple
import random

def parse_ion(name):
    ion=namedtuple('Ion', 'name anum mass charge')
    ion.name=name
    if name.upper()=='NA':
        ion.anum=11
        ion.mass=22.99
        ion.charge=1

    elif name.upper()=='K':
        ion.anum=19
        ion.mass=39.10
        ion.charge=1

    elif name.upper()=='CL':
        ion.anum=17
        ion.mass=35.45
        ion.charge=-1

    else:
        raise ValueError, "Unrecognized ion type '%s'" % name

    return ion

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


def Neutralize(mol, cation='NA', anion='CL', 
        chain='ION', chain2='ION2',
        solute_pad=5.0,
        ion_pad=3.0,
        concentration=0.0,
        verbose=False):

    """
    neutralize() -> replace water molecules with ions
  
    Optional arguments (default):
      cation  ('NA') -- species of cation.  Must be NA or K
      anion   ('CL') -- species of anion.  Currently only CL is supported
      prefix ('ION') -- counterions (those with the opposite charge of the
                          solute) are placed in segment ION; the other
                          species is put in segment ION2.

      solute_pad (5) -- minimum distance between placed ions and non-water.
      ion_pad    (3) -- minimum distance between placed ions.
      conc       (0) -- molar concentration of added ions with same charge
                          as the solute.  The number of such ions will be
                          given by 

                            nother = int((conc / 55.345) * (ntotalwat - nions))

                          where ntotalwat is the total number of waters in
                          the original system and nions is the number of 
                          counterions needed to achieve charge neutrality.

    """

    # first, we add sufficient counterions to neutralize.  Then we add 
    # counterions and counter-counterions until counter-counterions are
    # up to the desired concentration.  
    solute=mol.atomselect('not resname %s %s' % (cation, anion))
    cg = sum(a.charge for a in solute)
    print "neutralize: solute charge=%s" % cg

    nions = int(math.fabs(cg)+0.5)
    iontype = anion
    othertype = cation
    if cg < 0:
        iontype = cation
    othertype = anion

    iontype = parse_ion(iontype)
    othertype = parse_ion(othertype)

    # find the water residues
    water = mol.atomselect(
            'water and (not hydrogen) and (not within %f of (not water))' 
            % solute_pad)
    residues = list(set(a.residue for a in water))
    nwat = len(residues)
    print "waters available to be replaced by ions:", nwat

    # compute number of ions already present
    nions_prev = len(mol.atomselect('atomicnumber %d and not bonded' % iontype.anum))
    nother_prev = len(mol.atomselect('atomicnumber %d and not bonded' % othertype.anum))

    print "Starting with %d %s ions" % (nions_prev, iontype.name)
    print "Starting with %d %s ions" % (nother_prev, othertype.name)

    # convert molar concentration to number based on available waters.  The
    # molar concentration of water is about 55.345 mol/L.  Use all available
    # waters to calculate the number of ions to add.
    ntotalwat = len(set(a.residue for a in mol.atomselect('water')))
    print "Starting with %d water molecules" % ntotalwat
    nother = int((concentration / 55.345) * (ntotalwat - nions + nions_prev))
    nions += nother

    # subtract off the ions already present in solution
    nions -= nions_prev
    nother -= nother_prev
    if nions < 0 or nother < 0:
        raise RuntimeError, "Too many ions already in solution"

    if nwat < nions:
        raise RuntimeError, "Only %d waters found; not enough to neutralize" % nwat

    if nother > 0:
        print "Adding %d %s ions" % (nother, othertype.name)
    if nions > 0:
        print "Adding %d %s ions" % (nions, iontype.name)

    # Shuffle the residues
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
          rj = residues[j]
          try: pj = centers[rj.id]
          except KeyError:
              pj = centers[rj] = compute_center(rj)
          d2 = dist2(pi,pj)
          if d2 < ionpad2:
              try: 
                  del residues[j]
              except: 
                  raise RuntimeError, "Not enough waters or too large ion_pad."
          else:
            j += 1

    if nions > 0:
        ionchain = mol.addChain() 
        ionchain.name = chain 
        for i in range(nions):
            res=residues[i]

            ionres=ionchain.addResidue()
            ionres.num=i+1

            ionatm=ionres.addAtom()
            ionatm.pos=compute_center(res)
            ionatm.name=iontype.name
            ionatm.atomic_number = iontype.anum
            ionatm.charge = iontype.charge
            ionatm.mass = iontype.mass

            res.remove()

    if nother > 0:
        otherchain = mol.addChain()
        otherchain.name = chain2
        for i in range(nions, nions+nother):
            res=residues[i]

            ionres=ionchain.addResidue()
            ionres.num=i+1-nions

            ionatm=ionres.addAtom()
            ionatm.pos=compute_center(res)
            ionatm.name=othertype.name
            ionatm.atomic_number = othertype.anum
            ionatm.charge = othertype.charge
            ionatm.mass = othertype.mass

            res.remove()

    mol.reassignGids()

