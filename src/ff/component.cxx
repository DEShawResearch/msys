#include "../ff.hxx"

using namespace desres::msys;

void ff::build(Tuples& tuples, SystemPtr mol, IdList const& fragment) {

    tuples.atoms = fragment;
    MultiIdList& bonds = tuples.bonds;
    MultiIdList& angles = tuples.angles;
    MultiIdList& dihedrals = tuples.dihedrals;

    for (auto ai : fragment) {
        if (mol->atomFAST(ai).atomic_number==0) continue;
        auto ibonded = mol->bondedAtoms(ai);
        for (auto aj : ibonded) {
            if (mol->atomFAST(aj).atomic_number==0) continue;
            for (auto ak : ibonded) {
                if (mol->atomFAST(ak).atomic_number==0) continue;
                if (aj < ak) {
                    angles.emplace_back(IdList{aj,ai,ak});
                }
            }
            if (ai > aj) continue;
            bonds.emplace_back(IdList{ai,aj});
            auto jbonded = mol->bondedAtoms(aj);
            for (auto ah : ibonded) {
                if (mol->atomFAST(ah).atomic_number==0) continue;
                if (ah==aj) continue;
                for (auto ak : jbonded) {
                    if (mol->atomFAST(ak).atomic_number==0) continue;
                    if (ak==ai || ak==ah) continue;
                    if (ah < ak) {
                        dihedrals.emplace_back(IdList{ah,ai,aj,ak});
                    } else {
                        dihedrals.emplace_back(IdList{ak,aj,ai,ah});
                    }
                }
            }
        }
    }
}

