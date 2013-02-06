#include "aromatic_new.hxx"
#include "aromatic.hxx"
#include "bondFilters.hxx"
#include <queue>

using namespace desres::msys;

namespace {

    /* Get the ring system containing this bond */
    void get_ring_system(Id bond, SystemPtr sys, IdList& ring_system_atoms,
            IdList& ring_system_bonds, MultiIdList& rings) {

        ring_system_atoms.clear();
        ring_system_bonds.clear();
        rings.clear();

        /* Generate map from bond to list of potentially aromatic SSSR rings
         * containing that bond */
        annotation_t annot(sys);
        boost::shared_ptr<MultiIdList> SSSR = sys->allRelevantSSSR();
        MultiIdList bond_to_rings(sys->maxBondId());
        for (unsigned i = 0; i < SSSR->size(); ++i) {
            const IdList& ring = SSSR->at(i);
            bool possibly_aromatic = true;
            for (unsigned j = 0; j < ring.size(); ++j) {
                if (annot.hybridization(ring[j]) != 2)
                    possibly_aromatic = false;
            }
            /* Ring is potentially aromatic if all atoms are sp2 hybridized */
            if (!possibly_aromatic) continue;
            for (unsigned j = 0; j < ring.size(); ++j) {
                Id ring_bond = sys->findBond(ring[j], ring[(j+1)%ring.size()]);
                bond_to_rings[ring_bond].push_back(i);
            }
        }
        if (bond_to_rings[bond].size() == 0)
            /* Bond is not in potentially aromatic ring */
            return;

        /* Get the ring system containing this bond */
        std::set<Id> atom_set;
        std::set<Id> bond_set;
        std::set<unsigned> ring_set;
        std::queue<Id> unprocessed_bonds;
        bond_set.insert(bond);
        atom_set.insert(sys->bond(bond).i);
        atom_set.insert(sys->bond(bond).j);
        unprocessed_bonds.push(bond);
        while (unprocessed_bonds.size() > 0) {
            Id front = unprocessed_bonds.front();
            unprocessed_bonds.pop();
            /* Loop through all potentially aromatic SSSR rings containing
             * bond 'front' */
            for (unsigned i = 0; i < bond_to_rings[front].size(); ++i) {
                ring_set.insert(bond_to_rings[front][i]);
                const IdList& ring = SSSR->at(bond_to_rings[front][i]);
                for (unsigned j = 0; j < ring.size(); ++j) {
                    Id ring_bond = sys->findBond(ring[j],
                            ring[(j+1)%ring.size()]);
                    /* If ring bond is new, add to unprocessed bond queue */
                    if (bond_set.insert(ring_bond).second) {
                        unprocessed_bonds.push(ring_bond);
                        atom_set.insert(sys->bond(ring_bond).i);
                        atom_set.insert(sys->bond(ring_bond).j);
                    }
                }
            }
        }

        /* Convert sets to lists */
        ring_system_atoms.insert(ring_system_atoms.end(), atom_set.begin(),
                atom_set.end());
        ring_system_bonds.insert(ring_system_bonds.end(), bond_set.begin(),
                bond_set.end());
        BOOST_FOREACH(unsigned i, ring_set)
            rings.push_back(SSSR->at(i));
    }

    /* Check if a ring or ring system satisfies Huckel's rule */
    bool check_aromatic(const IdList& atoms, const IdList& bonds,
            SystemPtr sys, const std::vector<bool>& arom) {
        bondedVirtualsFilter filter(sys);
        std::set<Id> atom_set(atoms.begin(), atoms.end());
        int electron_count = 0;
        /* If double bond, add 2 to electron count */
        for (unsigned i = 0; i < bonds.size(); ++i)
            electron_count += 2 * (sys->bond(bonds[i]).order - 1);
        for (unsigned i = 0; i < atoms.size(); ++i) {
            IdList bonds_for_atom = sys->filteredBondsForAtom(atoms[i], filter);
            bool has_double = false;
            for (unsigned j = 0; j < bonds_for_atom.size(); ++j) {
                if (sys->bond(bonds_for_atom[j]).order == 2) {
                    has_double = true;
                    Id other = sys->bond(bonds_for_atom[j]).other(atoms[i]);
                    /* If external double bond to aromatic atom; add 1 to
                     * electron count */
                    if (arom.size() > 0 && arom[other] &&
                            atom_set.find(other) == atom_set.end())
                        ++electron_count;
                }
                else if (sys->bond(bonds_for_atom[j]).order != 1)
                    MSYS_FAIL("sp2 ring atom has bond order != 1 or 2");
            }
            /* We know atom is sp2; if no double bonds, must have lone pair in
             * ring--add 2 to electron count */
            if (!has_double)
                electron_count += 2;
        }
        /* Use Huckel's rule */
        return (electron_count % 4 == 2);
    }
}

bool desres::msys::IsAromaticBond(SystemPtr sys, Id bond) {
    IdList ring_system_atoms;
    IdList ring_system_bonds;
    MultiIdList rings;
    get_ring_system(bond, sys, ring_system_atoms, ring_system_bonds, rings);

    /* Check if entire ring system is aromatic */
    std::vector<bool> atom_aromatic(sys->maxAtomId(), false);
    if (check_aromatic(ring_system_atoms, ring_system_bonds, sys,
            atom_aromatic))
        return true;

    /* Check if individual rings are aromatic */
    std::vector<bool> ring_aromatic(rings.size(), false);
    bool detected = true;
    /* Do while the previous iteration marked at least one new aromatic ring */
    while (detected) {
        detected = false;
        for (unsigned i = 0; i < rings.size(); ++i) {
            /* If ring has already been marked as aromatic, leave as aromatic */
            if (ring_aromatic[i]) continue;
            IdList ring_bonds(rings[i].size());
            bool in_this_ring = false;
            for (unsigned j = 0; j < rings[i].size(); ++j) {
                ring_bonds[j] = sys->findBond(rings[i][j],
                        rings[i][(j+1)%rings[i].size()]);
                if (ring_bonds[j] == bond)
                    in_this_ring = true;
            }
            if (check_aromatic(rings[i], ring_bonds, sys, atom_aromatic)) {
                if (in_this_ring)
                    /* A ring containing this bond is aromatic */
                    return true;
                detected = true;
                ring_aromatic[i] = true;
                /* Update aromatic atoms for next iteration */
                for (unsigned j = 0; j < rings[i].size(); ++j)
                    atom_aromatic[rings[i][j]] = true;
            }
        }
    }
    /* All aromatic rings have been detected; none of them contain this bond */
    return false;
}

bool desres::msys::IsAromaticAtom(SystemPtr sys, Id atom) {
    bondedVirtualsFilter filter(sys);
    IdList bonds_for_atom = sys->filteredBondsForAtom(atom, filter);
    for (unsigned i = 0; i < bonds_for_atom.size(); ++i)
        if (IsAromaticBond(sys, bonds_for_atom[i]))
            return true;
    return false;
}
