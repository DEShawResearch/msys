

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

sys.setrecursionlimit(10000)

def det(t0,t1,int,norm):
    v0 = t1-t0
    v1 = int-t0
    c = numpy.cross(v0,v1)
    det = numpy.dot(c, norm)
    return det

def check_line_intersect_tri(t0,t1,t2,l0,l1):
    # print "t0=",t0,"t1=",t1,"t2=",t2
    l01 = l1-l0
    # print "l01=", l01

    t01 = t1-t0
    t12 = t2-t1

    norm = numpy.cross(t01,t12)
    # print "t01=",t01,"t12=",t12,"norm=", norm
    norm_dot_line = numpy.dot(l01,norm)
    
    if norm_dot_line == 0: 
        #print "line is parallel to triangle"
        return False

    l0mt0 = l0-t0
    # print "l0mt0",l0mt0
    t = -numpy.dot(norm,l0mt0) / numpy.dot(norm,l01)
    if t<=0 or t>=1:
        #print "intersection is past the line's endpoints"
        return False

    intersection = l0 + l01*t
    # print "intersection=",intersection
    
    signed_areas = numpy.array([det(t0,t1,intersection,norm),
                                det(t1,t2,intersection,norm),
                                det(t2,t0,intersection,norm)])
    # print "signed areas: ",signed_areas/2/sum(norm**2)**0.5
    cw = numpy.sign(signed_areas)
    # print "clockwise: ",cw
    return cw[0]!=0 and cw[0]==cw[1] and cw[0]==cw[2]

def ut_intersection():
    t0 = numpy.array([0.0, 0.0, 0.0]);
    t1 = numpy.array([0.0, 3.0, 3.0]);
    t2 = numpy.array([0.0, 2.0, 0.0]);

    l0 = numpy.array([-1.0, 1.0, 1.0]);
    l1n = numpy.array([1.0, 0.5, 1.0]);
    l1p = numpy.array([1.0, 1.5, 1.0]);
    
    assert check_line_intersect_tri(t0,t1,t2,l0,l1n) == False
    assert check_line_intersect_tri(t0,t1,t2,l0,l1p) == True

#######################################################################
#######################################################################
#######################################################################
#######################################################################

def FindKnots(mol, max_cycle_size=None, selection='all', verbose=False):

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
                if check_line_intersect_tri(
                        cp[0    ], cp[idx  ], cp[idx+1], pos[ai], pos[aj]):
                    num_int += 1
                    if verbose: print "==> intersection:",cycle,bond
                    results.append((cycle, bond))
    if verbose: print "Total intersections: %d" % num_int
    return results


