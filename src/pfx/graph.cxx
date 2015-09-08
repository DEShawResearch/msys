#include "graph.hxx"

using namespace desres::msys;
using namespace desres::msys::pfx;

Graph::Graph(SystemPtr mol) 
: g(mol->maxAtomId()), edgecount(0) {
    if (mol->maxAtomId() != mol->atomCount()) {
        MSYS_FAIL("System has deleted atoms, so graph would be incorrect.");
    }
    for (Id i=0, n=mol->maxAtomId(); i<n; i++) {
        IdList atoms = mol->bondedAtoms(i);
        for (Id j=0, m=atoms.size(); j<m; j++) {
            add_edge(i,atoms[j]);
        }
    }
}

