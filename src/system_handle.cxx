#include "system_handle.hxx"
#include "atomsel.hxx"
#include "append.hxx"
#include "clone.hxx"
#include "dms.hxx"
#include "mae.hxx"

#include <fstream>

using namespace desres::msys;


AtomList SystemHandle::atoms() {
    return Atom::list(_system, _system->atoms());
}

BondList SystemHandle::bonds() {
    return Bond::list(_system, _system->bonds());
}

ResidueList SystemHandle::residues() {
    return Residue::list(_system, _system->residues());
}

ChainList SystemHandle::chains() {
    return Chain::list(_system, _system->chains());
}

Atom SystemHandle::atom(Id id) {
    return Atom(_system, id);
}

Bond SystemHandle::bond(Id id) {
    return Bond(_system, id);
}

Residue SystemHandle::residue(Id id) {
    return Residue(_system, id);
}

Chain SystemHandle::chain(Id id) {
    return Chain(_system, id);
}

Atom SystemHandle::addAtom() {
    return Atom::create(_system);
}

Bond SystemHandle::addBond(Atom a1, Atom a2) {
    return Bond::create(_system, a1, a2);
}

Residue SystemHandle::addResidue() {
    return Residue::create(_system);
}

Chain SystemHandle::addChain() {
    return Chain::create(_system);
}

Bond SystemHandle::findBond(Atom a1, Atom a2) {
    Id id;
    if (_system->findBond(a1.id(), a2.id(), &id)) {
        return bond(id);
    }
    return bond(BadId);
}

SystemHandle desres::msys::CreateSystem() {
    return System::create();
}

SystemHandle 
desres::msys::LoadMAE(String const& path, bool ignore_unrecognized ) {
    return desres::msys::ImportMAE(path,ignore_unrecognized);
}

SystemHandle 
desres::msys::LoadDMS(String const& path, bool with_forcefield) {
    return desres::msys::ImportDMS(path,with_forcefield);
}

void desres::msys::SaveDMS(SystemHandle h, String const& path) {
    desres::msys::ExportDMS(h.ptr(), path);
}

void desres::msys::SaveMAE(SystemHandle h, String const& path) {
    desres::msys::ExportMAE(h.ptr(), path);
}

void SystemHandle::delAtoms(const AtomList& atoms) {
    _system->delAtoms(Atom::ids(atoms, _system));
}
void SystemHandle::delBonds(const BondList& bonds) {
    _system->delBonds(Bond::ids(bonds, _system));
}
void SystemHandle::delResidues(const ResidueList& residues) {
    _system->delResidues(Residue::ids(residues, _system));
}
void SystemHandle::delChains(const ChainList& chains) {
    _system->delChains(Chain::ids(chains, _system));
}

AtomList SystemHandle::atomselect( String const& seltext ) const {
    return Atom::list(ptr(), Atomselect(ptr(), seltext));
}

AtomList SystemHandle::appendSystem( SystemHandle src ) {
    return Atom::list(ptr(), AppendSystem(ptr(), src.ptr()));
}

SystemHandle desres::msys::CloneSystem( AtomList const& atoms ) {
    if (!atoms.size()) {
        return CreateSystem();
    }
    IdList ids;
    SystemPtr p = Atom::extract(atoms, ids);
    return desres::msys::Clone(p, ids);
}

