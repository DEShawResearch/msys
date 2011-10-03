#include "view.hxx"

using namespace desres::msys;

View::View(System* sys, const IdList& ids) 
: _sys(sys) {
    _atoms.insert(ids.begin(), ids.end());
    for (IdList::const_iterator i=ids.begin(), e=ids.end(); i!=e; ++i) {
        Id residue = _sys->atom(*i).residue;
        Id chain = _sys->residue(residue).chain;
        _residues.insert(residue);
        _chains.insert(chain);
        IdList bonds = _sys->bondsForAtom(*i);
        for (IdList::const_iterator b=bonds.begin(); b!=bonds.end(); ++b) {
            bond_t& bond = _sys->bond(*b);
            if (_atoms.count(bond.i) && _atoms.count(bond.j)) {
                _bonds.insert(*b);
            }
        }
    }
}

