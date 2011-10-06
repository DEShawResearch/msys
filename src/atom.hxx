#ifndef desres_msys_atom_hxx
#define desres_msys_atom_hxx

#include "ref.hxx"

namespace desres { namespace msys {

    class Atom;
    class Bond;
    class Residue;
    
    typedef std::vector<Bond> BondList;
    typedef std::vector<Atom> AtomList;
    
    class Atom : public Ref<Atom> {

        atom_t& atom() { return sys()->atom(id()); }
        const atom_t& atom() const { return sys()->atom(id()); }

    public:
        Atom(SystemPtr sys=SystemPtr(), Id id=BadId) : Ref<Atom>(sys,id) {}
        static Atom create( SystemPtr sys );

        ValueRef prop(const String& key) {
            return sys()->atomPropValue(id(), key);
        }
    
        Id bondCount() const { return sys()->bondCountForAtom(id()); }
        BondList bonds() const;

        AtomList bondedAtoms() const;
    
        Bond addBond( Atom a );

        /* accessors */
        Id gid() const { return atom().gid; }
        void setGid(Id gid) { atom().gid=gid; }

        Id fragid() const { return atom().fragid; }

        const Float* pos() const { return &atom().x; }
        void setPos(const Float* pos) {
            std::copy(pos,pos+3,&atom().x);
        }

        const Float* vel() const { return &atom().vx; }
        void setVel(const Float* vel) {
            std::copy(vel,vel+3,&atom().vx);
        }

        Float charge() const { return atom().charge; }
        void setCharge(Float charge) { atom().charge=charge; }

        Float mass() const { return atom().mass; }
        void setMass(Float mass) { atom().mass=mass; }

        int atomicNumber() const { return atom().atomic_number; }
        void setAtomicNumber(int anum) { atom().atomic_number=anum; }

        int formalCharge() const { return atom().formal_charge; }
        void setFormalCharge(int charge) { atom().formal_charge=charge; }

        String name() const { return atom().name; }
        void setName(const String& name) { atom().name=name; }

        Residue residue() const;
        void setResidue(Residue res);
    };

}}

#endif
