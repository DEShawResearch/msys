#include "builder.hxx"
#include "../../src/clone.hxx"
#include <boost/foreach.hpp>
#include <stdexcept>


namespace desres { namespace msys { namespace builder {

    namespace {
        Id find_atom(SystemPtr mol, Id res, Id prev, Id next, adef_t const& def) {
            if (fabs(def.rel)>1) {
                throw std::runtime_error(
                        "unsupported relative bond with magnitude > 1");
            }
            res = def.rel== 1 ? next :
                  def.rel==-1 ? prev :
                  res;
            if (!bad(res)) {
                BOOST_FOREACH(Id atm, mol->atomsForResidue(res)) {
                    if (def.name == mol->atom(atm).name) return atm;
                }
            }
            return BadId;
        }
    }
    void build( defs_t const& defs, SystemPtr mol, Id chain ) {
    
        mol->delBonds(mol->bonds());

        std::vector<resdef_t> resdefs;
        IdList const& residues = mol->residuesForChain(chain);
    
        /* assign resdef to each residue */
        BOOST_FOREACH(Id ires, residues) {
            residue_t& res = mol->residue(ires);
            ResdefMap::const_iterator idef = defs.resdefs.find(res.name);
            if (idef==defs.resdefs.end()) {
                fprintf(stderr, "No topology for residue '%s : %d'\n",
                        res.name.c_str(), res.resid);
                return;
            }
            resdef_t resdef = idef->second;
            /* TODO: apply patch to resdef if first or last */

            /* Remove unrepresented atoms */
            IdList atoms = mol->atomsForResidue(ires);
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
            BOOST_FOREACH(AtomMap::value_type const& adef, resdef.atoms) {
                Id atm = BadId;
                if (!amap.count(adef.first)) {
                    printf("added atom %s\n", adef.first.c_str());
                    atm = mol->addAtom(ires);
                } else {
                    atm = amap[adef.first];
                }
                mol->atom(atm).name = adef.first;
                mol->atom(atm).charge = adef.second.charge;
            }

            /* cache resdef for later */
            resdefs.push_back(resdef);
        }

        /* add bonds.  We to have added all atoms before we can do this,
         * hence the need for a second pass. */
        for (Id i=0; i<residues.size(); i++) {
            resdef_t const& resdef = resdefs[i];
            Id res=residues[i];
            Id prev = i==0 ? BadId : residues[i-1];
            Id next = i==residues.size()-1 ? BadId : residues[i+1];

            BOOST_FOREACH(bond_t const& b, resdef.bonds) {
                Id atm1 = find_atom(mol, res, prev, next, b.def1);
                Id atm2 = find_atom(mol, res, prev, next, b.def2);
                if (bad(atm1) || bad(atm2)) {
                    printf("Skipping bond to nonexistent atom %s-%s\n",
                            b.def1.name.c_str(), b.def2.name.c_str());
                    continue;
                }
                mol->addBond(atm1, atm2);
            }
        }

        /* clone and replace */
        SystemPtr tmp = Clone(mol, mol->atoms());
        mol.swap(tmp);
    }
}}}
