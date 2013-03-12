#include "annotated_system.hxx"
#include "elements.hxx"
#include "sssr.hxx"
#include <queue>
#include <assert.h>

desres::msys::AnnotatedSystem::AnnotatedSystem(SystemPtr sys)
: _sys(sys), _atoms(sys->maxAtomId()), _bonds(sys->maxBondId()) {

    /* Assign valence, degree, hcount, lone_pairs, and preliminary
     * hybridization */
    BOOST_FOREACH(Id ai, sys->atoms()) {
        atom_data_t& a = _atoms[ai];
        double res_valence = 0;
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
            res_valence += bnd.resonant_order;
            a.degree += 1;
        }
        int formal_charge = sys->atom(ai).formal_charge;
        int val = DataForElement(sys->atom(ai).atomic_number).nValence;
        int electrons = val - a.valence - formal_charge;
        if (electrons < 0 || electrons % 2)
            MSYS_FAIL("Invalid formal charge or bond orders for atom "
                    << ai << " of system " << sys->name);
        a.lone_pairs = electrons / 2;
        a.resonant_lone_pairs = (val - res_valence
                - sys->atom(ai).resonant_charge) / 2;
        /* Check that resonant lone pairs is not negative */
        if (a.resonant_lone_pairs < -0.00001)
            MSYS_FAIL("Invalid resonant charge or bond orders for atom "
                    << ai << " of system " << sys->name);
        if (sys->atom(ai).atomic_number == 1 || a.degree == 0)
            a.hybridization = 0;
        else
            /* Add 0.00001 before taking floor to prevent possible floating
             * point error */
            a.hybridization = std::min(3, std::max(1,
                    int(a.degree + a.resonant_lone_pairs + 0.00001 - 1)));
    }
    /* Get rings and ring systems, assign ring_bonds and rings_idx */
    compute_ring_systems();
    /* Assign aromatic */
    compute_aromaticity();
    /* If sp3 boron/carbon/nitrogen/oxygen AND have lone pairs AND
     * (aromatic OR bonded to atom with double bonds): become sp2 */
    BOOST_FOREACH(Id ai, sys->atoms()) {
        atom_data_t& a = _atoms[ai];
        int element = sys->atom(ai).atomic_number;
        /* Check against 0.00001 instead of 0 to prevent possible floating
         * point error */
        if (element >= 5 && element <= 8 && a.hybridization == 3
                && a.resonant_lone_pairs > 0.00001) {
            if (a.aromatic) {
                a.hybridization = 2;
                continue;
            }
            BOOST_FOREACH(Id aj, sys->bondedAtoms(ai)) {
                if (_atoms[aj].hybridization == 2)
                    a.hybridization = 2;
            }
        }
    }
}

void desres::msys::AnnotatedSystem::compute_ring_systems() {
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

    MultiIdList SSSR_possibly_aromatic;
    IdList rid_map;
    for (unsigned i = 0; i < SSSR.size(); ++i) {
        bool possibly_aromatic = true;
        /* All atoms in aromatic ring must be potentially sp2 (meaning sp2 or
         * hybridization greater than 2 with free electrons) */
        BOOST_FOREACH(Id atom, SSSR[i]) {
            /* Check against 0.00001 instead of 0 to prevent possible floating
             * point error */
            if (_atoms[atom].hybridization != 2 &&
                    (_atoms[atom].hybridization < 2
                     || _atoms[atom].resonant_lone_pairs < 0.00001))
                possibly_aromatic = false;
        }
        if (possibly_aromatic) {
            SSSR_possibly_aromatic.push_back(SSSR[i]);
            rid_map.push_back(i);
        }
    }

    MultiIdList ring_systems = RingSystems(_sys, SSSR_possibly_aromatic);
    BOOST_FOREACH(const IdList& ring_sys, ring_systems) {
        std::set<Id> atom_set;
        std::set<Id> bond_set;
        std::set<Id> ring_set;
        BOOST_FOREACH(Id rid, ring_sys) {
            const IdList& ring = SSSR_possibly_aromatic[rid];
            for (unsigned i = 0; i < ring.size(); ++i) {
                atom_set.insert(ring[i]);
                Id bond = _sys->findBond(ring[i], ring[(i+1)%ring.size()]);
                bond_set.insert(bond);
            }
            ring_set.insert(rid_map[rid]);
        }
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
        /* If no double bonds and has lone pair, lone pair is in ring--add 2 to
         * electron count */
        if (!has_double && _atoms[atom].lone_pairs > 0) electron_count += 2;
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
