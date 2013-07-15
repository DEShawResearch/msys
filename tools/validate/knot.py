

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

def ut_intersection():
    t0 = numpy.array([0.0, 0.0, 0.0]);
    t1 = numpy.array([0.0, 3.0, 3.0]);
    t2 = numpy.array([0.0, 2.0, 0.0]);

    l0 = numpy.array([-1.0, 1.0, 1.0]);
    l1n = numpy.array([1.0, 0.5, 1.0]);
    l1p = numpy.array([1.0, 1.5, 1.0]);
    
    assert msys.LineIntersectsTriangle(l0,l1n,t0,t1,t2) == False
    assert msys.LineIntersectsTriangle(l0,l1p,t0,t1,t2) == True

ut_intersection()

#######################################################################
#######################################################################
#######################################################################
#######################################################################

def FindKnots(mol, max_cycle_size=None, selection='all', verbose=False):

    from msys import LineIntersectsTriangle
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

    # fetch positions
    pos = mol.getPositions()
    ids = set(a.id for a in atoms)
    bonds=[]
    for b in mol.bonds:
        ai, aj = b.first.id, b.second.id
        if ai in ids and aj in ids:
            bonds.append((ai,aj))

    results = []
    num_int = 0

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
                if LineIntersectsTriangle(
                        pos[ai], pos[aj], cp[0], cp[idx], cp[idx+1]):
                    num_int += 1
                    if verbose: print "==> intersection:",cycle,bond
                    results.append((cycle, bond))
    if verbose: print "Total intersections: %d" % num_int
    return results


