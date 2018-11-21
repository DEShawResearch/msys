#include "pfx/pfx.hxx"
#include <assert.h>
#include <stdio.h>

using namespace desres::msys::pfx;

static const unsigned bonds[][2] = {
    {0,1}, {0,2},
    {4,3}, {4,6}, {5,6},
    {7,8}, {9,8}, {9,10}, {10,11}, {11,7}
};
static const unsigned nbonds = sizeof(bonds)/sizeof(bonds[0]);

static float pos[12][3] = {
    // first fragment
    {0.0,2,0},
    {0.0,4.1,0},
    {0.1, 2.2, .3},
    // second fragment
    {3,1,3},
    {3,1.1,3},
    {3.1, 1.2, 3.3},
    {3.2, 1.0, 6.5},
    // third fragment
    {1,0,2},
    {1,0.8,2},
    {1.1, 4.2, 2.3},
    {0.8, 0.0, 2.5},
    {0.9, 0.0, 2.1}
};
static const unsigned natoms = sizeof(pos)/sizeof(pos[0]);

static float box[] = {
    0.70710678, 0.70710678, 0.,
    -1.41421356, 1.41421356, 0.,
    0, 0, 3
};

static const unsigned glue[] = {1,8,9};
static const unsigned nglue = sizeof(glue)/sizeof(glue[0]);

int main() {
    Graph g(natoms);
    assert(g.nverts()==natoms);
    assert(g.nedges()==0);
    for (unsigned i=0; i<nbonds; i++) {
        assert(g.add_edge(bonds[i][0], bonds[i][1]));
    }
    assert(g.nverts()==natoms);
    assert(g.nedges()==nbonds);

    std::vector<Graph::Id> sorted_bonds;
    g.copy_bonds(std::back_inserter(sorted_bonds));
    assert(sorted_bonds.size()==2*nbonds);

    std::vector<Graph::Id> comps, sizes;
    g.copy_components(std::back_inserter(comps),
                      std::back_inserter(sizes));
    assert(sizes.size()==3);
    assert(comps.size()==natoms);
    assert(sizes[0]==3);
    assert(sizes[1]==4);
    assert(sizes[2]==5);

    Pfx pfx(g, true);
    pfx.glue(nglue, glue);
    pfx.apply(&pos[0][0], box, (float *)NULL);
    //for (unsigned i=0; i<natoms; i++) {
        //printf("%2u: %8.5f %8.5f %8.5f\n", i,
                //pos[i][0], pos[i][1], pos[i][2]);
    //}

    return 0;
}
