
import os, msys, math
from collections import namedtuple

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

def Neutralize(mol, cation='NA', anion='CL', 
        chain='ION', counterchain='ION2',
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
    residues = set(a.residue for a in water)
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


    if nwat < nions:
        raise RuntimeError, "Only %d waters found; not enough to neutralize" % nwat

    if nother > 0:
        print "Adding %d %s ions" % (nother, othertype.name)
    if nions > 0:
        print "Adding %d %s ions" % (nions, iontype.name)

    return

    # Import the structure into psfgen
    s = psfgen.structure()
    s.import_structure(psf)
    s.import_coordinates(pdb)

    # Shuffle the residues
    residues = numpy.random.permutation(residues).tolist()
    mass = atomsel('residue %d' % residues[0]).get('mass')

    # Remove overly close residues among the first nions + nother waters
    # Cache the centers to save time
    centers={}
    ionpad2 = ion_pad * ion_pad
    for i in range(nions + nother):
        ri = residues[i]
        try: pi = centers[ri]
        except KeyError:
            pi = centers[ri] = numpy.array(atomsel('residue %d' % ri).center(mass))
        j = i+1
        while j < nions+nother:
          rj = residues[j]
          try: pj = centers[rj]
          except KeyError:
              pj = centers[rj] = numpy.array(atomsel('residue %d' % rj).center(mass))

          d2 = sum((pi-pj)**2)
          if d2 < ionpad2:
              try: del residues[j]
              except: raise RuntimeError, "Not enough waters or too large ion_pad."
          else:
            j += 1

    if nions > 0:
      # Delete the first nions of them, after fetching their positions
      pos = vmdnumpy.positions()
      positions = []
      for i in range(nions):
        wat = atomsel('residue %d' % residues[i])
        segid = wat.get('segid')[0]
        resid = wat.get('resid')[0]
        avgpos = wat.center(wat.get('mass'))

        # store the position and delete the water
        positions.append(avgpos)
        s.delatom(segid, resid)

      # Create a new segment for the counterions
      s.build( prefix, residues = [(resid, iontype) for resid in range(nions)] )

      # Set positions for the ions
      for resid in range(nions):
          s.set_coordinates( prefix, resid, iontype, list(positions[resid]) )

    # Create a segment for the counter-counterions
    if nother > 0:
        positions = []
        for i in range(nions, nions+nother):
          wat = atomsel('residue %d' % residues[i])
          segid = wat.get('segid')[0]
          resid = wat.get('resid')[0]
          avgpos = wat.center(wat.get('mass'))

          # store the position and delete the water
          positions.append(avgpos)
          s.delatom(segid, resid)

        prefix2 = '%s2' % prefix
        s.build(prefix2, residues=[(resid, othertype) for resid in range(nother)])
        for resid in range(nother):
          s.set_coordinates(prefix2, resid, othertype, list(positions[resid]) )

    if nions < 0 or nother < 0:
        # delete ions and replace them with waters.  Don't bother shuffling the
      # list of ion positions.
      pos = vmdnumpy.positions()
      sel=atomsel('resname %s' % iontype)
      if len(sel) < -nions:
          raise RuntimeError, "Cannot decrease concentration - not enough ions to delete!"
      ion_inds=sel.get('index')[:-nions]
      ion_segs=sel.get('segid')[:-nions]
      ion_resids=sel.get('resid')[:-nions]
      sel=atomsel('resname %s' % othertype)
      if len(sel) < -nother:
          raise RuntimeError, "Cannot decrease concentration - not enough ions to delete!"
      other_inds=sel.get('index')[:-nother]
      other_segs=sel.get('segid')[:-nother]
      other_resids=sel.get('resid')[:-nother]
      inds=ion_inds + other_inds
      segs=ion_segs + other_segs
      resids=ion_resids + other_resids

      print "Replacing %d %s ions and %d %s ions with TIP3" % (
              -nions, iontype, -nother, othertype)

      positions=vmdnumpy.positions()[inds]
      for segid, resid in zip(segs, resids):
          s.delatom(segid, resid)

      # create a segment for the new waters
      prefix3 = '%s3' % prefix
      s.build(prefix3, residues=[(resid, 'TIP3') for resid in range(len(inds))])
      for resid in range(len(inds)):
          s.set_coordinates(prefix3, resid, 'OH2', list(positions[resid]) )
      s.guess_coordinates()

    # Write out the new version
    s.writepsf('%s.psf' % output, cmap=False)
    s.writepdb('%s.pdb' % output )

    print "neutralize finished normally."

