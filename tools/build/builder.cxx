#include "builder.hxx"
#include <boost/foreach.hpp>

namespace desres { namespace msys { namespace builder {

    void build( defs_t const& defs, SystemPtr mol, Id chain ) {
    
        std::vector<resdef_t> resdefs;
        IdList const& residues = mol->residuesForChain(chain);
    
        /* assign resdef to each residue */
        for (Id i=0; i<residues.size(); i++) {
            residue_t& res = mol->residue(residues[i]);
            ResdefMap::const_iterator idef = defs.resdefs.find(res.name);
            if (idef==defs.resdefs.end()) {
                fprintf(stderr, "No topology for residue '%s : %d'\n",
                        res.name.c_str(), res.resid);
                return;
            }
            resdef_t resdef = idef->second;

            /* TODO: apply patch to resdef if first or last */

            /* Remove unrepresented atoms */
            IdList atoms = mol->atomsForResidue(residues[i]);
            std::map<std::string, Id> amap;
            for (Id j=0; j<atoms.size(); j++) {
                msys::atom_t& atm = mol->atom(atoms[j]);
                AtomMap::const_iterator adef=resdef.atoms.find(atm.name);
                /* remove if not in resdef or if duplicate name */
                if (adef==resdef.atoms.end() || amap.count(atm.name)) {
                    printf("deleted atom %s\n", atm.name.c_str());
                    mol->delAtom(atoms[j]);
                } else {
                    amap[atm.name] = atoms[j];
                }
            }
            /* Add missing atoms */
            BOOST_FOREACH(AtomMap::value_type adef, resdef.atoms) {
                Id atm = BadId;
                if (!amap.count(adef.first)) {
                    printf("added atom %s\n", adef.first.c_str());
                    atm = mol->addAtom(residues[i]);
                } else {
                    atm = amap[adef.first];
                }
                mol->atom(atm).name = adef.first;
                mol->atom(atm).charge = adef.second.charge;
            }

            /* cache resdef for later */
            resdefs.push_back(resdef);
        }

    }
}}}
