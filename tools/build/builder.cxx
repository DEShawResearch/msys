#include "builder.hxx"

namespace desres { namespace msys { namespace builder {

    void build( defs_t const& defs, SystemPtr mol, Id chain ) {
    
        /* assign resdef to each residue */
        std::vector<resdef_t> resdefs;
        IdList const& residues = mol->residuesForChain(chain);
    
        for (Id i=0; i<residues.size(); i++) {
            residue_t& res = mol->residue(residues[i]);
            ResdefMap::const_iterator idef = defs.resdefs.find(res.name);
            if (idef==defs.resdefs.end()) {
                fprintf(stderr, "No topology for residue '%s : %d'\n",
                        res.name.c_str(), res.resid);
                continue;
            }
            resdefs.push_back(idef->second);
        }
    
    }
}}}
