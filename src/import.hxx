#ifndef desres_msys_import_hxx
#define desres_msys_import_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* A helper object for importers that manages the logic of 
     * mapping atoms to residues and residue to chains */
    class SystemImporter {
        SystemPtr sys;

        struct ChnKey {
            Id     ct;
            String name;
            String segid;

            ChnKey() {}
            ChnKey(Id c, String const& nm, String const& seg)
            : ct(c), name(nm), segid(seg) {}

            bool operator<(ChnKey const& c) const {
                if (ct!=c.ct) return ct<c.ct;
                int rc = name.compare(c.name);
                if (rc) return rc<0;
                return segid.compare(c.segid)<0;
            }
        };
                
        struct ResKey {
            Id      chain;
            int     resnum;
            String  resname;
            String  insertion;

            ResKey() {}
            ResKey(Id chn, int num, String const& name, String const& insert) 
            : chain(chn), resnum(num), resname(name), insertion(insert) {}

            bool operator<(const ResKey& r) const {
                if (chain!=r.chain) return chain<r.chain;
                if (resnum!=r.resnum) return resnum<r.resnum;
                if (resname!=r.resname) return resname<r.resname;
                return insertion.compare(r.insertion)<0;
            }
        };

        typedef std::map<ChnKey,Id> ChnMap;
        typedef std::map<ResKey,Id> ResMap;
        ResMap resmap;
        ChnMap chnmap;

        Id chnid;
        Id resid;

    public:
        explicit SystemImporter(SystemPtr s) 
        : sys(s), chnid(BadId), resid(BadId) {}

        /* process existing atoms in the system */
        void initialize(IdList const& atoms);

        /* mark a chain as terminated, so that subsequent atoms with
         * the same chain name will be added to a new chain object. */
        bool terminateChain(std::string chain, std::string segid, Id ct=0);

        /* add an atom, after first constructing necessary parent
         * chain and/or residue object.  All string inputs will
         * have leading and trailing whitespace removed. */
        Id addAtom(std::string chain, std::string segid, 
                   int resnum, std::string resname, 
                   std::string atomname,
                   std::string insertion="",
                   Id ct=0);
    };

}}

#endif
