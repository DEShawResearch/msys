#include "annotated_system.hxx"
#include "elements.hxx"
#include "sssr.hxx"
#include <queue>

using namespace desres::msys;

AnnotatedSystem::AnnotatedSystem(SystemPtr sys, unsigned flags)
: _atoms(sys->maxAtomId()), _bonds(sys->maxBondId()) {

    for (Id b : sys->bonds()) {
        bond_t const& bnd = sys->bondFAST(b);
        Id i = bnd.i;
        Id j = bnd.j;
        int ni = sys->atomFAST(i).atomic_number;
        int nj = sys->atomFAST(j).atomic_number;
        if (ni<1 || nj<1) continue;
        if (bnd.order<1) {
            MSYS_FAIL("Invalid bond order for bond " << b << " atoms "
                    << i << "," << j << " of system " << sys->name);
        }
        atom_data_t& ai = _atoms[i];
        atom_data_t& aj = _atoms[j];
        if (ai.degree==atom_data_t::MaxBonds) {
            MSYS_FAIL("Too many bonds (" << atom_data_t::MaxBonds
                    << ") for atom " << i << " of system " << sys->name);
        }
        if (aj.degree==atom_data_t::MaxBonds) {
            MSYS_FAIL("Too many bonds (" << atom_data_t::MaxBonds
                    << ") for atom " << j << " of system " << sys->name);
        }
        ai.valence += bnd.order;
        aj.valence += bnd.order;
        ai.bond[ai.degree++] = b;
        aj.bond[aj.degree++] = b;
        if (ni==1) aj.hcount++;
        if (nj==1) ai.hcount++;

        bond_data_t& bij = _bonds[b];
        bij.i = i;
        bij.j = j;
        bij.order = bnd.order;
    }

    int nrad=0;
    for (Id i : sys->atoms()) {
        atom_t const& atm = sys->atomFAST(i);
        int anum = atm.atomic_number;
        if (anum<1) continue;

        atom_data_t& a = _atoms[i];

        int formal_charge = atm.formal_charge;
        int val = DataForElement(anum).nValence;
        int aval = DataForElement(anum).additionalValence;
        int electrons = val - a.valence - formal_charge;
        if (electrons < 0 && aval > 0) {
            electrons = aval - a.valence - formal_charge;
        }
        if (electrons < 0) {
            std::stringstream ss;
            ss << "Invalid formal charge or bond orders for atom "
                << AbbreviationForElement(anum) << " "
                << i << " of system " << sys->name;
            if (flags & AllowBadCharges) {
                _errors.push_back(ss.str());
            } else MSYS_FAIL(ss.str());
        }
        nrad+=electrons%2;
        a.lone_electrons=electrons;
        a.formal_charge = formal_charge;
        a.atomic_number = anum;

        if (anum==1 || a.degree==0) {
            a.hybridization = 0;
        } else {
            a.hybridization = std::max(1, a.degree+(a.lone_electrons+1)/2 - 1);
        }
    }
    if (nrad >1) {
        std::stringstream ss;
        ss << "Invalid formal charge or bond orders ( "<<nrad <<
                      " radical centers detected ) for system " << sys->name;
        if (flags & AllowBadCharges) {
            _errors.push_back(ss.str());
        } else MSYS_FAIL(ss.str());
    }

    /* Get rings and ring systems, assign ring_bonds and rings_idx */
    compute_ring_systems(sys);
    /* Assign aromatic */
    compute_aromaticity(sys);
    /* If sp3 AND have lone electrons AND group={14,15,16} AND
     * (aromatic OR bonded to atom with double bonds): become sp2 */
    for (Id ai : sys->atoms()) {
        atom_data_t& a = _atoms[ai];
        int group=GroupForElement(sys->atom(ai).atomic_number);
        bool validGroup=group>=14 && group<=16;
        if (a.hybridization == 3 && a.lone_electrons > 1 && validGroup ) {
            if (a.aromatic) {
                a.hybridization = 2;
                continue;
            }
            for (Id aj : sys->bondedAtoms(ai)){
                if ( _atoms[aj].hybridization == 2){
                    a.hybridization = 2;
                    break;
                }
                if(a.degree==1 && _atoms[aj].hybridization == 3){
                    for (Id bi : sys->bondsForAtom(aj)){
                        bond_t const& bnd = sys->bond(bi);
                        if(bnd.order != 2) continue;
                        if(_atoms[bnd.other(aj)].degree==1){
                            a.hybridization = 2;
                            break;
                        }
                    }
                    if(a.hybridization == 2) break;
                }
            }
        }
    }
}

IdList AnnotatedSystem::atoms() const {
    IdList ids;
    for (Id i=0, n=_atoms.size(); i<n; i++) {
        if (_atoms[i].atomic_number>0) ids.push_back(i);
    }
    return ids;
}

void AnnotatedSystem::compute_ring_systems(SystemPtr _sys) {
    MultiIdList SSSR = GetSSSR(_sys, _sys->atoms(), true);
    std::vector<IdSet> ring_bonds(_sys->maxAtomId());
    for (const IdList& ring : SSSR) {
        _rings.push_back(ring_t());
        _rings.back().atoms = ring;
        for (unsigned i = 0; i < ring.size(); ++i) {
            _atoms.at(ring[i]).rings_idx.push_back(_rings.size()-1);
            Id bond = _sys->findBond(ring[i], ring[(i+1)%ring.size()]);
            _rings.back().bonds.push_back(bond);
            _bonds.at(bond).rings_idx.push_back(_rings.size()-1);
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
         * sp3 with free electrons) */
        for (Id atom : SSSR[i]) {
            if (! ( ( _atoms[atom].hybridization == 3 && _atoms[atom].lone_electrons > 1) ||
                    _atoms[atom].hybridization == 2))
                possibly_aromatic = false;
        }
        if (possibly_aromatic) {
            SSSR_possibly_aromatic.push_back(SSSR[i]);
            rid_map.push_back(i);
        }
    }

    MultiIdList ring_systems = RingSystems(_sys, SSSR_possibly_aromatic);
    for (const IdList& ring_sys : ring_systems) {
        std::set<Id> atom_set;
        std::set<Id> bond_set;
        std::set<Id> ring_set;
        for (Id rid : ring_sys) {
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

bool AnnotatedSystem::is_aromatic(SystemPtr _sys, const IdList& atoms, const IdList& bonds) {

    int electron_count = 0;
    for (Id bond : bonds) {
        /* Internal double bond, add 2 to electron count */
        if (_sys->bond(bond).order == 2)
            electron_count += 2;
    }
    std::set<Id> atom_set(atoms.begin(), atoms.end());
    for (Id atom : atoms) {
        IdList bondedAtoms = _sys->bondedAtoms(atom);
        bool has_double = false;
        for (Id bonded : bondedAtoms) {
            int bondedAtomicNumber=_sys->atom(bonded).atomic_number;
            if (bondedAtomicNumber < 1) continue;
            Id bond = _sys->findBond(atom, bonded);
            if (_sys->bond(bond).order == 2) {
                has_double = true;
                if (atom_set.find(bonded)== atom_set.end()){
                    if(_atoms[bonded].aromatic){
                       /* External double bond to aromatic atom; add 1 to
                        * electron count */
                        ++electron_count;
                    }else if(_sys->atom(atom).atomic_number==6 && bondedAtomicNumber==6){
                       /* External C=C bond to nonaromatic carbon destroys aromaticity */
                        return false;
                    }
                }
            }
        }
        /* If no double bonds and has lone pair, lone pair is in ring--add 2 to
         * electron count */
        if (!has_double && _atoms[atom].lone_electrons > 1) electron_count += 2;
    }
    /* Use Huckel's rule */
    return (electron_count % 4 == 2);
}

void AnnotatedSystem::compute_aromaticity(SystemPtr _sys) {
    bool detected = true;
    /* Do while previous iteration marked a new aromatic atom or bond */
    while (detected) {
        std::vector<bool> ring_aromatic(_rings.size(), false);
        for (const ring_system_t& ring_sys : _ring_systems) {
            /* Check if entire ring system is aromatic */
            if (is_aromatic(_sys, ring_sys.atoms, ring_sys.bonds)) {
                for (Id ring : ring_sys.rings)
                    ring_aromatic[ring] = true;
            } else {
                /* Check if individual rings are aromatic */
                for (Id ring : ring_sys.rings) {
                    if (is_aromatic(_sys, _rings[ring].atoms, _rings[ring].bonds))
                        ring_aromatic[ring] = true;
                }
            }
        }
        detected = false;
        for (unsigned i = 0; i < _rings.size(); ++i) {
            if (ring_aromatic[i]) {
                /* Set aromaticity for atoms/bonds of this ring */
                for (Id atom : _rings[i].atoms) {
                    detected |= (!_atoms[atom].aromatic);
                    _atoms[atom].aromatic = true;
                }
                for (Id bond : _rings[i].bonds) {
                    detected |= (!_bonds[bond].aromatic);
                    _bonds[bond].aromatic = true;
                }
            } else { /* We do not undo previously detected aromaticity */ }
        }
    }
}
