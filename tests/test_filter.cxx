#include "load.hxx"
#include "atomsel.hxx"
#include <stdio.h>

using namespace desres::msys;

namespace {
    struct GetH {
        SystemPtr _m;
        explicit GetH(SystemPtr m) : _m(m) {}
        bool operator()(atom_t const& a) const {
            return a.atomic_number==1;
        }
        bool operator()(bond_t const& b) const {
            return _m->atom(b.i).atomic_number==1 ||
                   _m->atom(b.j).atomic_number==1;
        }
    };

    struct GetHBond {
        SystemPtr _m;
        Id ai;
        GetHBond(SystemPtr m, Id a) : _m(m), ai(a) {}
        bool operator()(bond_t const& b) const {
            return _m->atom(b.other(ai)).atomic_number==1;
        }
    };
}

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        SystemPtr mol = Load(argv[i]);
        IdList ids = Atomselect(mol, "protein and noh");
        for (Id j=0; j<ids.size(); j++) {
            Id ai = ids[j];
            IdList hbonds = mol->filteredBondsForAtom(ai, GetH(mol));
            IdList hatoms = mol->filteredBondedAtoms(ai, GetH(mol));
            printf("atom %u %s has %u bonds, %lu of which are to H\n",
                    ai, mol->atom(ai).name.c_str(), mol->bondCountForAtom(ai),
                    hbonds.size());
            for (Id k=0; k<hatoms.size(); k++) {
                Id aj = hatoms[k];
                printf("  %u %s\n", aj, mol->atom(aj).name.c_str());
            }
            assert(mol->filteredBondsForAtom(ai,GetHBond(mol,ai))==hbonds);
        }
    }
    return 0;
}

