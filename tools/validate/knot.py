

"""
Given a system, find any instances where a bond is 'threaded' through a ring, 
like so

                /-------\
       \       /         \
        \     /           \
         \----\  ---------------\
               \         /       \
                \-------/         \

The algorithm works as follows:

    1. Produce a list of all cycles in the bond topology (i.e. rings)
    2. For each ring:
        a. Use boxing and distance cutoffs to reduce the number of bonds to check against
        b. Divide the ring into N triangles
        c. Check for a triangle-line intersection between the triangle and each relevant bond
"""
# Copyright (C) 2010 D.E. Shaw Research
# @author Adam Lerer
# @author Justin Gullingsrud (converted to msys)


import numpy
import sys
import math
import msys
from msys.wrap import Wrapper

from time import time

def fooLineIntersectsTriangle(*args): return False

def ut_intersection():
    t0 = numpy.array([0.0, 0.0, 0.0]);
    t1 = numpy.array([0.0, 3.0, 3.0]);
    t2 = numpy.array([0.0, 2.0, 0.0]);

    l0 = numpy.array([-1.0, 1.0, 1.0]);
    l1n = numpy.array([1.0, 0.5, 1.0]);
    l1p = numpy.array([1.0, 1.5, 1.0]);
    
    assert msys.LineIntersectsTriangle(l0,l1n,t0,t1,t2) == False
    assert msys.LineIntersectsTriangle(l0,l1p,t0,t1,t2) == True

#######################################################################
#######################################################################
#######################################################################
#######################################################################

def FindKnots(mol, max_cycle_size=None, selection='all', verbose=False):

    from msys import LineIntersectsTriangle
    t0=time()

    mol = mol.clone()
    if selection is None: selection = 'all'
    atoms = mol.select(selection) 
    cycles = msys.GetSSSR(atoms)
    if verbose: print "Found %d cycles" % len(cycles)
    if not cycles: return []
    if verbose: print "Largest cycle has size %d" % max(map(len,cycles))
    cycles = [tuple(a.id for a in c) for c in cycles]
    if max_cycle_size is not None:
        max_cycle_size = int(max_cycle_size)
        cycles = [c for c in cycles if len(c) <= max_cycle_size]
        if verbose:
            print "Reduced to %d cycles of length <= %d" % (
                    len(cycles), max_cycle_size)
  
    t1=time()

    # fetch positions
    ids = set(a.id for a in atoms)
    bonds=[]
    for b in mol.bonds:
        ai, aj = b.first.id, b.second.id
        if ai in ids and aj in ids:
            bonds.append((ai,aj))

    t2=time()

    results = list()
    found = set()
    box = mol.getCell()
    wrapper = Wrapper(mol)

    t3=time()

    # apply shift of half the box size in all directions to catch knots that
    # cross periodic boundaries.
    for ishift in range(2):
        if ishift>0:
            mol.translate(0.5*(box[0]+box[1]+box[2]))
        wrapper.wrap()
        pos = mol.getPositions()

        for icycle, cycle in enumerate(cycles):
            cycle_sel='index ' + ' '.join(map(str,cycle))
            sel = '(%s) and exwithin 10 of %s' % (selection, cycle_sel)
            ids = set(mol.selectIds(sel))
            cp = [pos[i] for i in cycle]
            cycle_inds = range(1,len(cycle)-1)
            for bond in bonds:
                ai, aj = bond
                if ai not in ids: continue
                for idx in cycle_inds:
                    key = (ai,aj,0,idx,idx+1)
                    if key in found: continue
                    if LineIntersectsTriangle(
                            pos[ai], pos[aj], cp[0], cp[idx], cp[idx+1]):
                        if verbose: print "==> intersection:",cycle,bond
                        results.append((cycle, bond))
                        found.add(key)
    if verbose: print "Total intersections: %d" % len(found)

    t4=time()

    if False:
        print "Time: %8.3f %8.3f %8.3f %8.3f" % (
                (t1-t0)*1000,
                (t2-t1)*1000,
                (t3-t2)*1000,
                (t4-t3)*1000,)

    return results


