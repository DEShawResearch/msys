#ifndef desres_msys_chain_hxx
#define desres_msys_chain_hxx

#include "ref.hxx"

namespace desres { namespace msys {

    class Residue;
    
    class Chain : public Ref<Chain> {

        chain_t& chain() { return sys()->chain(id()); }
        const chain_t& chain() const { return sys()->chain(id()); }

    public:
        Chain(SystemPtr sys=SystemPtr(), Id id=BadId) : Ref<Chain>(sys,id) {}
        static Chain create( SystemPtr ent );
        Id size() const { return sys()->residueCountForChain(id()); }

        /* accessors */
        String name() const { return chain().name; }
        void setName(const String& name) { chain().name=name; }
    
        Residue addResidue();
        ResidueList residues() const {
            return Residue::list(sys(), sys()->residuesForChain(id()));
        }
    };
    
    typedef std::vector<Chain> ChainList;
}}

#endif
