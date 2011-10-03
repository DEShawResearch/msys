#ifndef desres_msys_bond_hxx
#define desres_msys_bond_hxx

#include "ref.hxx"

namespace desres { namespace msys {

    class Atom;
    
    class Bond : public Ref<Bond> {

        bond_t& bond() { return sys()->bond(id()); }
        const bond_t& bond() const { return sys()->bond(id()); }

    public:
        Bond(SystemPtr sys=SystemPtr(), Id id=BadId) : Ref<Bond>(sys,id) {}
        static Bond create( SystemPtr sys, Atom a1, Atom a2 );

        /* accessors */
        Float order() const { return bond().order; }
        void setOrder(Float order) { bond().order=order; }
        ValueRef prop(const String& key) {
            Id col = sys()->bondProps()->propIndex(key);
            return sys()->bondProps()->value(id(), col);
        }

        Id i() const { return bond().i; }
        Id j() const { return bond().j; }
    
        Atom first() const { return Atom(sys(), bond().i); }
        Atom second() const { return Atom(sys(), bond().j); }
        Atom other(Atom atom) const { 
            return Atom(sys(), bond().other(atom.id()));
        }

        AtomList atoms() const;
    };
    
}}

#endif
