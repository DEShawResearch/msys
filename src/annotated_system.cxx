#include "annotated_system.hxx"
#include "elements.hxx"
#include "sssr.hxx"
#include <queue>
#include <assert.h>

desres::msys::AnnotatedSystem::AnnotatedSystem(SystemPtr sys, unsigned flags)
: _sys(sys), _atoms(sys->maxAtomId()), _bonds(sys->maxBondId()) {

    /* Assign valence, degree, hcount, lone_electrons, and preliminary
     * hybridization */
    int nrad=0;
    for (Id ai : sys->atoms()) {
        int anum = sys->atomFAST(ai).atomic_number;
        if (anum < 1) continue;
        atom_data_t& a = _atoms[ai];
        for(Id bi : sys->bondsForAtom(ai)) {
            bond_t const& bnd = sys->bond(bi);
            Id aj = bnd.other(ai);
            int anum_j = sys->atom(aj).atomic_number;
            if (anum_j < 1) continue;
            if (bnd.order == 0)
                MSYS_FAIL("Invalid bond order for bond "
                        << bi << " atoms " << ai << "," << aj << " of system " << sys->name);
            if (anum_j == 1) ++a.hcount;
            a.valence += bnd.order;
            a.degree += 1;
            _bonds[bi].order = bnd.order;
        }
        int formal_charge = sys->atom(ai).formal_charge;
        int val = DataForElement(anum).nValence;
        int aval = DataForElement(anum).additionalValence;
        int electrons = val - a.valence - formal_charge;
        if (electrons < 0 && aval > 0) {
            electrons = aval - a.valence - formal_charge;
        }
        if (electrons < 0) {
            //fprintf(stderr, "Atom %s nval %d q %d val %d\n", 
                    //AbbreviationForElement(anum), val, formal_charge,a.valence);
            try {
                MSYS_FAIL("Invalid formal charge or bond orders for atom "
                    << AbbreviationForElement(anum) << " "
                    << ai << " of system " << sys->name);
            }
            catch (std::exception &e) {
                if (flags & AllowBadCharges) {
                    _errors.push_back(e.what());
                } else {
                    throw;
                }
            }
        }
        nrad+=electrons%2;
        a.lone_electrons=electrons;
        a.formal_charge = formal_charge;
        a.atomic_number = anum;

        if (anum == 1 || a.degree == 0){
            a.hybridization = 0;
        }else{
            a.hybridization = std::max(1, a.degree + (a.lone_electrons+1)/2 - 1);
        }
    }
    if (nrad >1) {
        try {
            MSYS_FAIL("Invalid formal charge or bond orders ( "<<nrad <<
                      " radical centers detected ) for system " << sys->name);
        }
        catch (std::exception& e) {
            if (flags & AllowBadCharges) {
                _errors.push_back(e.what());
            } else {
                throw;
            }
        }
    }

    /* Get rings and ring systems, assign ring_bonds and rings_idx */
    compute_ring_systems();
    /* Assign aromatic */
    compute_aromaticity();
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
                    for (Id bi : _sys->bondsForAtom(aj)){
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
         * sp3 with free electrons) */
        BOOST_FOREACH(Id atom, SSSR[i]) {
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

bool desres::msys::AnnotatedSystem::is_aromatic(const IdList& atoms,
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

void desres::msys::AnnotatedSystem::compute_aromaticity() {
    bool detected = true;
    /* Do while previous iteration marked a new aromatic atom or bond */
    while (detected) {
        std::vector<bool> ring_aromatic(_rings.size(), false);
        BOOST_FOREACH(const ring_system_t& ring_sys, _ring_systems) {
            /* Check if entire ring system is aromatic */
            if (is_aromatic(ring_sys.atoms, ring_sys.bonds)) {
                BOOST_FOREACH(Id ring, ring_sys.rings)
                    ring_aromatic[ring] = true;
            } else {
                /* Check if individual rings are aromatic */
                BOOST_FOREACH(Id ring, ring_sys.rings) {
                    if (is_aromatic(_rings[ring].atoms, _rings[ring].bonds))
                        ring_aromatic[ring] = true;
                }
            }
        }
        detected = false;
        for (unsigned i = 0; i < _rings.size(); ++i) {
            if (ring_aromatic[i]) {
                /* Set aromaticity for atoms/bonds of this ring */
                BOOST_FOREACH(Id atom, _rings[i].atoms) {
                    detected |= (!_atoms[atom].aromatic);
                    _atoms[atom].aromatic = true;
                }
                BOOST_FOREACH(Id bond, _rings[i].bonds) {
                    detected |= (!_bonds[bond].aromatic);
                    _bonds[bond].aromatic = true;
                }
            } else { /* We do not undo previously detected aromaticity */ }
        }
    }
}
