#include "annotated_system.hxx"
#include "elements.hxx"
#include "sssr.hxx"
#include <queue>
#include <assert.h>

desres::msys::AnnotatedSystem::AnnotatedSystem(SystemPtr sys)
: _sys(sys), _atoms(sys->maxAtomId()), _bonds(sys->maxBondId()) {

    /* Assign valence, degree, hcount, lone_pairs */
    BOOST_FOREACH(Id ai, sys->atoms()) {
        atom_data_t& a = _atoms[ai];
        BOOST_FOREACH(Id bi, sys->bondsForAtom(ai)) {
            bond_t const& bnd = sys->bond(bi);
            Id aj = bnd.other(ai);
            int anum_j = sys->atom(aj).atomic_number;
            if (anum_j < 1) continue;
            if (bnd.order == 0)
                MSYS_FAIL("Invalid bond order for bond "
                        << bi << " of system " << sys->name);
            if (anum_j == 1) ++a.hcount;
            a.valence += bnd.order;
            a.degree += 1;
        }
        int formal_charge = sys->atom(ai).formal_charge;
        int electrons = DataForElement(sys->atom(ai).atomic_number).nValence
            - a.valence - formal_charge;
        if (electrons < 0 || electrons % 2)
            MSYS_FAIL("Invalid formal charge or bond orders for atom "
                    << ai << " of system " << sys->name);
        a.lone_pairs = electrons / 2;
    }
    /* Get rings and ring systems, assign ring_bonds and rings_idx */
    compute_SSSR_rings();
    compute_ring_systems();
    /* Assign aromatic */
    compute_aromaticity();
    /* Assign hybridization */
    BOOST_FOREACH(Id ai, sys->atoms()) {
        atom_data_t& a = _atoms[ai];
        a.hybridization = a.degree + a.lone_pairs - 1;
        /* If sp3 AND have lone pairs AND (aromatic OR are bonded to C or N
         * that has double bonds): become sp2 */
        if (a.hybridization == 3 && a.lone_pairs > 0) {
            if (a.aromatic) {
                a.hybridization = 2;
                continue;
            }
            BOOST_FOREACH(Id aj, sys->bondedAtoms(ai)) {
                int nj = sys->atom(aj).atomic_number;
                if (!(nj==6 || nj==7)) continue;
                BOOST_FOREACH(Id bj, sys->bondsForAtom(aj))
                    if (sys->bond(bj).order==2) a.hybridization = 2;
            }
        }
    }
}

void desres::msys::AnnotatedSystem::compute_SSSR_rings() {
    MultiIdList SSSR = GetSSSR(_sys, _sys->atoms(), true);
    std::vector<IdSet> ring_bonds(_sys->maxAtomId());
    BOOST_FOREACH(const IdList& ring, SSSR) {
        _rings.push_back(ring_t());
        _rings.back().atoms = ring;
        for (unsigned i = 0; i < ring.size(); ++i) {
            assert(ring[i] < _atoms.size());
            _atoms[ring[i]].rings_idx.push_back(_rings.size()-1);
            Id bond = _sys->findBond(ring[i], ring[(i+1)%ring.size()]);
            _rings.back().bonds.push_back(bond);
            assert(bond < _bonds.size());
            _bonds[bond].rings_idx.push_back(_rings.size()-1);
            ring_bonds[ring[i]].insert(bond);
            ring_bonds[ring[(i+1)%ring.size()]].insert(bond);
        }
    }
    for (unsigned i = 0; i < ring_bonds.size(); ++i)
        if (ring_bonds[i].size() > 0)
            _atoms[i].ring_bonds = ring_bonds[i].size();
}

void desres::msys::AnnotatedSystem::compute_ring_systems() {
    /* Map from bond to list of potentially aromatic SSSR rings containing that
     * bond */
    MultiIdList bond_to_rings(_sys->maxBondId());
    for (unsigned i = 0; i < _rings.size(); ++i) {
        bool possibly_aromatic = true;
        BOOST_FOREACH(Id atom, _rings[i].atoms) {
            /* Atom must be sp2, or have hybridization > 2 with a lone pair */
            int hyb = _atoms[atom].degree + _atoms[atom].lone_pairs - 1;
            if (hyb != 2 && (hyb < 2 || _atoms[atom].lone_pairs == 0))
                possibly_aromatic = false;
        }
        if (!possibly_aromatic) continue;
        BOOST_FOREACH(Id bond, _rings[i].bonds)
            bond_to_rings[bond].push_back(i);
    }

    std::vector<bool> processed_bonds(_sys->maxBondId(), false);
    BOOST_FOREACH(Id bond, _sys->bonds()) {
        if (processed_bonds[bond]) continue;
        processed_bonds[bond] = true;
        if (bond_to_rings[bond].size() == 0) continue;
        /* Get the ring system containing this bond */
        std::set<Id> atom_set;
        std::set<Id> bond_set;
        std::set<Id> ring_set;
        std::queue<Id> unprocessed_bonds;
        bond_set.insert(bond);
        atom_set.insert(_sys->bond(bond).i);
        atom_set.insert(_sys->bond(bond).j);
        unprocessed_bonds.push(bond);
        while (unprocessed_bonds.size() > 0) {
            Id front = unprocessed_bonds.front();
            unprocessed_bonds.pop();
            /* Loop through all potentially aromatic SSSR rings containing
             * bond 'front' */
            BOOST_FOREACH(Id ring, bond_to_rings[front]) {
                ring_set.insert(ring);
                BOOST_FOREACH(Id ring_bond, _rings[ring].bonds) {
                    /* If ring bond is new, add to unprocessed bond queue */
                    if (bond_set.insert(ring_bond).second) {
                        unprocessed_bonds.push(ring_bond);
                        atom_set.insert(_sys->bond(ring_bond).i);
                        atom_set.insert(_sys->bond(ring_bond).j);
                        processed_bonds[ring_bond] = true;
                    }
                }
            }
        }
        /* Add this ring system */
        _ring_systems.push_back(ring_system_t());
        _ring_systems.back().atoms = IdList(atom_set.begin(), atom_set.end());
        _ring_systems.back().bonds = IdList(bond_set.begin(), bond_set.end());
        _ring_systems.back().rings = IdList(ring_set.begin(), ring_set.end());
    }
}

bool desres::msys::AnnotatedSystem::flag_aromatic(const IdList& atoms,
        const IdList& bonds) {
    int electron_count = 0;
    BOOST_FOREACH(Id bond, bonds) {
        /* Internal double bond, add 2 to electron count */
        if (_sys->bond(bond).order == 2)
            electron_count += 2;
    }
    std::set<Id> atom_set(atoms.begin(), atoms.end());
    BOOST_FOREACH(Id atom, atoms) {
        IdList bondedAtoms = _sys->bondedAtoms(atom);
        bool has_double = false;
        BOOST_FOREACH(Id bonded, bondedAtoms) {
            if (_sys->atom(bonded).atomic_number < 1) continue;
            Id bond = _sys->findBond(atom, bonded);
            if (_sys->bond(bond).order == 2) {
                has_double = true;
                if (_atoms[bonded].aromatic && atom_set.find(bonded)
                        == atom_set.end())
                    /* External double bond to aromatic atom; add 1 to
                     * electron count */
                    ++electron_count;
            }
        }
        /* If no double bonds, must have lone pair in ring--add 2 to
         * electron count */
        if (!has_double) electron_count += 2;
    }
    /* Use Huckel's rule */
    if (electron_count % 4 == 2) {
        BOOST_FOREACH(Id atom, atoms)
            _atoms[atom].aromatic = true;
        BOOST_FOREACH(Id bond, bonds)
            _bonds[bond].aromatic = true;
        return true;
    }
    return false;
}

void desres::msys::AnnotatedSystem::compute_aromaticity() {
    BOOST_FOREACH(const ring_system_t& ring_sys, _ring_systems) {
        /* Check if entire ring system is aromatic */
        if (flag_aromatic(ring_sys.atoms, ring_sys.bonds))
            continue;
        /* Check if individual rings are aromatic */
        std::vector<bool> ring_aromatic(ring_sys.rings.size(), false);
        bool detected = true;
        /* Do while the previous iteration marked at least one new aromatic ring */
        while (detected) {
            detected = false;
            for (unsigned i = 0; i < ring_sys.rings.size(); ++i) {
                /* If ring has already been marked as aromatic, leave as aromatic */
                if (ring_aromatic[i]) continue;
                if (flag_aromatic(_rings[ring_sys.rings[i]].atoms,
                            _rings[ring_sys.rings[i]].bonds)) {
                    detected = true;
                    ring_aromatic[i] = true;
                }
            }
        }
    }
}
