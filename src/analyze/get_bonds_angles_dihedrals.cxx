#include "../analyze.hxx"

namespace desres { namespace msys {

void GetBondsAnglesDihedrals(SystemPtr sys,
        const IdList& atoms, std::vector<IdList>& non_pseudo_bonds,
        std::vector<IdList>& pseudo_bonds, std::vector<IdList>& angles,
        std::vector<IdList>& dihedrals) {

    /* To check that there are no bonds to external atoms */
    std::vector<bool> in_frag(sys->maxAtomId(), false);
    for (unsigned i = 0; i < atoms.size(); ++i)
        in_frag[atoms[i]] = true;

    non_pseudo_bonds.clear();
    pseudo_bonds.clear();
    angles.clear();
    dihedrals.clear();

    IdList bond(2);
    IdList angle(3);
    IdList dihedral(4);
    for (unsigned i = 0; i < atoms.size(); ++i) {
        Id ai = atoms[i];
        if (sys->atom(ai).atomic_number == 0)
            continue;
        const IdList& ibonded = sys->bondedAtoms(ai);
        for (unsigned j = 0, m = ibonded.size(); j < m; ++j) {
            Id aj = ibonded[j];
            if (!in_frag[aj])
                MSYS_FAIL("Cannot get tuples: incomplete fragment");
            if (sys->atom(aj).atomic_number == 0) {
                /* Add pseudo bond ai-aj */
                bond[0] = ai;
                bond[1] = aj;
                pseudo_bonds.push_back(bond);
                continue;
            }
            /* Add angles with center ai */
            for (unsigned k = 0; k < m; ++k) {
                Id ak = ibonded[k];
                if (!in_frag[ak])
                    MSYS_FAIL("Cannot get tuples: incomplete fragment");
                if (sys->atom(ak).atomic_number != 0 && aj < ak) {
                    angle[0] = aj;
                    angle[1] = ai;
                    angle[2] = ak;
                    angles.push_back(angle);
                }
            }
            if (ai > aj) continue;
            /* Add non-pseudo bond ai-aj */
            bond[0] = ai;
            bond[1] = aj;
            non_pseudo_bonds.push_back(bond);
            /* Add dihedrals with center ai-aj */
            const IdList& jbonded = sys->bondedAtoms(aj);
            for (unsigned h = 0; h < m; ++h) {
                Id ah = ibonded[h];
                if (!in_frag[ah])
                    MSYS_FAIL("Cannot get tuples: incomplete fragment");
                if (sys->atom(ah).atomic_number == 0 || ah == aj)
                    continue;
                for (unsigned k = 0; k < jbonded.size(); ++k) {
                    Id ak = jbonded[k];
                    if (!in_frag[ak])
                        MSYS_FAIL("Cannot get tuples: incomplete fragment");
                    if (sys->atom(ak).atomic_number == 0
                            || ak == ai || ak == ah)
                        continue;
                    if (ah < ak) {
                        dihedral[0] = ah;
                        dihedral[1] = ai;
                        dihedral[2] = aj;
                        dihedral[3] = ak;
                    }
                    else {
                        dihedral[0] = ak;
                        dihedral[1] = aj;
                        dihedral[2] = ai;
                        dihedral[3] = ah;
                    }
                    dihedrals.push_back(dihedral);
                }
            }
        }
    }
}
}}
