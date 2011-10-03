#include "atom.hxx"
#include "bond.hxx"
#include "residue.hxx"
#include "chain.hxx"

using namespace desres::msys;

namespace desres { namespace msys {
    template <> void Ref<Atom>::destroy()    { sys()->delAtom(id()); }
    template <> void Ref<Bond>::destroy()    { sys()->delBond(id()); }
    template <> void Ref<Residue>::destroy() { sys()->delResidue(id()); }
    template <> void Ref<Chain>::destroy()   { sys()->delChain(id()); }

    template <> bool Ref<Atom>::exists() const { 
        return sys()->hasBond(id()); 
    }
    template <> bool Ref<Bond>::exists() const { 
        return sys()->hasBond(id()); 
    }
    template <> bool Ref<Residue>::exists() const { 
        return sys()->hasResidue(id()); 
    }
    template <> bool Ref<Chain>::exists() const { 
        return sys()->hasChain(id()); 
    }

    template<> std::string Ref<Atom>::objtype() const { return "Atom"; }
    template<> std::string Ref<Bond>::objtype() const { return "Bond"; }
    template<> std::string Ref<Residue>::objtype() const { return "Residue"; }
    template<> std::string Ref<Chain>::objtype() const { return "Chain"; }

}}
    
Atom Atom::create(SystemPtr cmp) {
    return Residue::create(cmp).addAtom();
}

Bond Bond::create(SystemPtr ent, Atom a1, Atom a2) {
    return Bond(ent, ent->addBond(a1.id(), a2.id()));
}

Residue Residue::create(SystemPtr ent) {
    return Chain::create(ent).addResidue();
}

Chain Chain::create(SystemPtr ent) {
    return Chain(ent, ent->addChain());
}

Bond Atom::addBond( Atom a ) {
    return Bond::create(sys(), *this, a);
}

BondList Atom::bonds() const {
    return Bond::list(sys(), sys()->bondsForAtom(id()));
}

AtomList Atom::bondedAtoms() const {
    return Atom::list(sys(), sys()->bondedAtoms(id()));
}

AtomList Bond::atoms() const {
    IdList ids(2);
    ids[0]=i();
    ids[1]=j();
    return Atom::list(sys(), ids);
}

Atom Residue::addAtom() {
    return Atom(sys(), sys()->addAtom(id()));
}

Residue Chain::addResidue() {
    return Residue(sys(), sys()->addResidue(id()));
}

Residue Atom::residue() const {
    return Residue(sys(), atom().residue);
}
void Atom::setResidue(Residue res) {
    if (sys()!=res.sys()) {
        throw std::runtime_error("Cannot assign Atom to Residue from a different System");
    }
    sys()->setResidue(id(), res.id());
}

Chain Residue::chain() const {
    return Chain(sys(), residue().chain);
}
