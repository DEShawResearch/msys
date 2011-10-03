#ifndef desres_msys_residue_hxx
#define desres_msys_residue_hxx

#include "ref.hxx"

namespace desres { namespace msys {

    class Atom;
    class Chain;
    typedef std::vector<Atom> AtomList;
    
    class Residue : public Ref<Residue> {

        residue_t& residue() { return sys()->residue(id()); }
        const residue_t& residue() const { return sys()->residue(id()); }

    public:
        Residue(SystemPtr sys=SystemPtr(), Id id=BadId) 
        : Ref<Residue>(sys,id) {}

        static Residue create( SystemPtr ent );

        /* Accessors. */
        String name() const { return residue().name; }
        void setName(const String& name) { residue().name=name; }

        int num() const { return residue().num; }
        void setNum(int num) { residue().num=num; }

        Chain chain() const;

        Id size() const { return sys()->atomCountForResidue(id()); }
        AtomList atoms() const {
            return Atom::list(sys(), sys()->atomsForResidue(id()));
        }
        Atom addAtom();

    };
    
    typedef std::vector<Residue> ResidueList;
}}

#endif
