

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

"""
find_all_cycles contributed by Mathias Laurin <Mathias Laurin AT gmail com>
"""

def find_cycle_to_ancestor(spanning_tree, node, ancestor):
    """
    Find a cycle containing both node and ancestor.
    """
    path = []
    while (node != ancestor):
        if node is None:
            return []
        path.append(node)
        node = spanning_tree[node]
    path.append(node)
    path.reverse()
    return path

def find_all_cycles(graph):
    """
    Find all cycles in the given graph.

    This function will return a list of lists of nodes, which form cycles in the
    graph or an empty list if no cycle exists.
    """

    def dfs(node):
        """
        Depth-first search subfunction.
        """
        visited.add(node)
        # Explore recursively the connected component
        for each in graph[node]:
            if each not in visited:
                spanning_tree[each] = node
                dfs(each)
            else:
                if (spanning_tree[node] != each):
                    cycle = find_cycle_to_ancestor(spanning_tree, node, each)
                    if cycle:
                        cycles.append(cycle)

    visited = set()         # List for marking visited and non-visited nodes
    spanning_tree = {}      # Spanning tree
    cycles = []

    # Algorithm outer-loop
    for each in graph:
        # Select a non-visited node
        if each not in visited:
            spanning_tree[each] = None
            # Explore node's connected component
            dfs(each)

    return cycles

#######################################################################
#######################################################################
#######################################################################
#######################################################################

def FindKnots(mol, max_cycle_size=None, use_boxing=True, include_h2o=False,
        verbose=False):

    bonds = {}
    bond_list = []
    
    h2o_particle_list = set()
    if not include_h2o:
        h2o_particle_list.update(a.id for a in mol.select('water'))
    
    for b in mol.bonds:
        a0, a1 = b.atoms
        p0, p1 = a0.id, a1.id
        if p0 in h2o_particle_list: continue
        if p0 not in bonds: bonds[p0] = []
        if p1 not in bonds: bonds[p1] = []
        bonds[p0] += [p1]
        bonds[p1] += [p0]
        bond_list += [(p0,p1)]
    cycles = find_all_cycles(bonds)
    cycles = [tuple(cycle) for cycle in cycles]
    if verbose: print "Found %d cycles" % len(cycles)
    if max_cycle_size:
        cycles = [cycle for cycle in cycles if len(cycle) <= max_cycle_size]
        if verbose:
            print "Reduced to %d cycles of length <= %d" % (len(cycles), max_cycle_size)

    
    # fetch positions
    pos = mol.getPositions()

    # l = (c.execute("select x from global_cell where id=1").fetchone()[0],
    #      c.execute("select y from global_cell where id=2").fetchone()[0],
    #      c.execute("select z from global_cell where id=3").fetchone()[0])
    # print l
    def insert(b,bi,boxes):
        atom_pos = pos[b[bi]]
        #idx = [math.floor(atom_pos[i] / l[i] * 2) for i in range(3)]
        idx = [1 if atom_pos[i] > 0 else 0  for i in range(3)]
        full_idx = 4*idx[0]+2*idx[1]+idx[2]
        # print b,bi,atom_pos,"in",idx
        boxes[full_idx].add(b)

    # WARNING: there are (literally) 'corner cases' that get missed by this boxing strategy
    bond_boxes = []
    for i in range(8): bond_boxes += [set()]
    cycle_boxes = []
    for i in range(8): cycle_boxes += [set()]
    if use_boxing:
        for b in bond_list:
            insert(b,0,bond_boxes)
            insert(b,1,bond_boxes)
        cycle_boxes = []
        for i in range(8): cycle_boxes += [set()]
        for c in cycles:
            for i in range(len(c)):
                insert(tuple(c),i, cycle_boxes)
    else:
        bond_boxes[0]  = set(bond_list)
        cycle_boxes[0] = set(cycles)

    if verbose:
        print "Total bonds", len(bond_list), "; bond boxes", [len(x) for x in bond_boxes]
        print "Total cycles", len(cycles), "; cycle boxes", [len(x) for x in cycle_boxes]

    results = []
    num_int = 0
    for box in range(8):
        if verbose:
            print "Checking box %d of %d" % (box+1, 8 if use_boxing else 1)
        for cycle in cycle_boxes[box]:
            cycle_set = set(cycle)
            #print "checking cycle", cycle
            for bond in bond_boxes[box]:
                bond_rep  = pos[bond[0]]
                cycle_rep = pos[cycle[0]]
                if sum((bond_rep-cycle_rep)**2) > 100: continue # skip if reps are more than 10A apart
                if bond[0] in cycle_set or bond[1] in cycle_set: continue

                cp = [pos[atom] for atom in cycle]
                bp = [pos[atom] for atom in bond]
                for cycle_idx in range(1,len(cycle)-1):
                    #print bond,"vs",(cycle[0],cycle[cycle_idx],cycle[cycle_idx+1])
                    intersect = check_line_intersect_tri(cp[0          ],
                                                         cp[cycle_idx  ],
                                                         cp[cycle_idx+1],
                                                         bp[0],
                                                         bp[1]
                                        )
                    if intersect:
                        num_int += 1
                        if verbose: print "==> intersection:",cycle,bond
                        results.append((cycle, bond))
                #print c1,c2
    if verbose: print "Total intersections: %d" % num_int
    return results


