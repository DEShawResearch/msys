#include "molecule.hxx"
#include "types.hxx"

using namespace desres::msys;

Molecule::Molecule(unsigned natoms, unsigned nbonds)
: _natoms(natoms), _nbonds(nbonds) {
    if (_natoms != (unsigned short)(natoms)) {
        MSYS_FAIL("natoms " << natoms << "too large");
    }
    if (_nbonds != (unsigned short)(nbonds)) {
        MSYS_FAIL("nbonds " << nbonds << "too large");
    }
    _atoms.reset(new Atom[natoms]);
    _bonds.reset(new Bond[nbonds]);
}

Molecule::Molecule(Molecule const& c)
: _name(c._name), _natoms(c._natoms), _nbonds(c._nbonds) {
    _atoms.reset(new Atom[_natoms]);
    _bonds.reset(new Bond[_nbonds]);
    memcpy(_atoms.get(), c._atoms.get(), _natoms*sizeof(Atom));
    memcpy(_bonds.get(), c._bonds.get(), _nbonds*sizeof(Bond));
}

Molecule& Molecule::operator=(Molecule const& c) {
    Molecule tmp(c);
    _atoms.swap(tmp._atoms);
    _bonds.swap(tmp._bonds);
    _natoms = tmp._natoms;
    _nbonds = tmp._nbonds;
    return *this;
}
