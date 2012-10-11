#include "graph.hxx"
#include "system.hxx"
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>

using namespace desres::msys;

/* Randomly permute ID list */
void permute(IdList& ids) {
    for (unsigned i = 0; i < ids.size(); ++i) {
        unsigned j = rand() % ids.size();
        Id tmp = ids[i];
        ids[i] = ids[j];
        ids[j] = tmp;
    }
}

void print(const IdList& ids) {
    fprintf(stdout, "[");
    for (unsigned i = 0; i < ids.size(); ++i)
        fprintf(stdout, "%d ", ids[i]);
    fprintf(stdout, "]");
}

void print(const std::vector<std::pair<Id, Id> >& ids) {
    fprintf(stdout, "[");
    for (unsigned i = 0; i < ids.size(); ++i)
        fprintf(stdout, "(%d %d) ", ids[i].first, ids[i].second);
    fprintf(stdout, "]");
}

int main(int argc,char **argv){
    unsigned natoms=5000, nexternals=300;
    if (argc==3) {
        natoms = atoi(argv[1]);
        nexternals = atoi(argv[2]);
    } else if (argc!=1) {
        fprintf(stderr, "Usage: %s [num_atoms num_externals]\n", argv[0]);
        return 1;
    }

    unsigned seed = time(NULL);
    //fprintf(stdout, "Seed: %d\n", seed);
    srand(seed);
    SystemPtr sysv[2] = {System::create(), System::create()};
    for (int sysid = 0; sysid < 2; ++sysid) {
        SystemPtr sys = sysv[sysid];
        Id chain = sys->addChain();
        Id residue = sys->addResidue(chain);
        for (unsigned i = 0; i < natoms; ++i) {
            /* Add random atoms with atomic number 1, 2, 3, 4, or 5 */
            Id atom = sys->addAtom(residue);
            sys->atom(atom).atomic_number = rand() % 5 + 1;
            int numatoms = sys->atomCount();
            if (i == 0)
                continue;
            /* For each atom, add between 1 and 4 bonds to random existing
             * atoms */
            Id other;
            do {
                other = rand() % numatoms;
            } while (other == atom);
            sys->addBond(atom, other);
            for (unsigned j = 0; j < 3; ++j) {
                other = rand() % numatoms;
                if (other != atom)
                    sys->addBond(atom, other);
            }
        }
        for (unsigned i = 0; i < nexternals; ++i) {
            /* Add 2 external atoms */
            Id atom = sys->addAtom(residue);
            sys->atom(atom).atomic_number = 0;
            int numatoms = sys->atomCount();
            Id other;
            do {
                other = rand() % numatoms;
            } while (other == atom);
            sys->addBond(atom, other);
        }
    }
    IdList atoms0 = sysv[0]->atoms();
    IdList atoms1 = sysv[1]->atoms();
    GraphPtr G0=Graph::create(sysv[0], atoms0);
    GraphPtr G1=Graph::create(sysv[1], atoms1);
    permute(atoms0);
    permute(atoms1);
    GraphPtr G0perm=Graph::create(sysv[0], atoms0);
    GraphPtr G1perm=Graph::create(sysv[1], atoms1);
    std::vector<std::pair<Id, Id> > output;
    assert(G0->match(G1, output) == false);
    assert(G0->match(G0perm, output) == true);
    for (unsigned i = 0; i < natoms; ++i)
        assert(output[i].first == output[i].second);
    assert(G1->match(G1perm, output) == true);
    for (unsigned i = 0; i < natoms; ++i)
        assert(output[i].first == output[i].second);
    printf("OK\n");

    return 0;
}
