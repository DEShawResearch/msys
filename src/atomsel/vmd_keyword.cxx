#include "vmd_keyword.hxx"
#include <stdexcept>

using namespace desres::msys::atomsel;

namespace desres { namespace msys { namespace atomsel { namespace {

    typedef std::map<Id,int> ResidueMap;
    void find_connected_residues( SystemPtr sys,
                                  ResidueMap& rmap,
                                  Id res,
                                  int id ) {
        rmap[res] = id;
        Id chn = sys->residue(res).chain;
        static const char * SG = "SG"; /* disulfide name */
        /* find residues connected to res */
        IdList atoms = sys->atomsForResidue(res);
        for (Id i=0; i<atoms.size(); i++) {
            Id atm = atoms[i];
            int atm_is_sg = sys->atom(atm).name==SG;
            IdList bonds = sys->bondsForAtom(atm);
            for (Id j=0; j<bonds.size(); j++) {
                Id anbr = sys->bond(bonds[j]).other(atm);
                const std::string& aname = sys->atom(anbr).name;
                Id nbr = sys->atom(anbr).residue;
                if ( nbr != res &&
                        !rmap.count(nbr) &&
                        sys->residue(nbr).chain==chn &&
                        !(atm_is_sg && aname==SG)) {
                    find_connected_residues( sys, rmap, nbr, id );
                }
            }
        }
    }

    class Fragment : public Keyword {
        SystemPtr _sys;
        ResidueMap rmap;
    public:
        Fragment(SystemPtr sys)  
        : Keyword("fragment", KEY_INT), _sys(sys) {
            /* in keeping with the VMD way of doing things, rather than just 
             * the connected subgraphs of the bonded atoms, this is a little 
             * more chemistry-specific.  We instead find subgraphs of 
             * connected residues within the same chain, and don't follow 
             * disulfide bridges.  */
            int id=0;
            IdList residues = _sys->residues();
            for (Id i=0; i<residues.size(); i++) {
                Id res = residues[i];
                if (rmap.count(res)) continue; /* already marked */
                find_connected_residues( _sys, rmap, res, id++ );
            }
        }
        void iget(const Selection& s, std::vector<Int>& v) const {
            for (Id i=0; i<s.size(); i++) {
                if (s[i]) {
                    Id res = _sys->atom(i).residue;
                    ResidueMap::const_iterator it=rmap.find(res);
                    if (it==rmap.end()) {
                        throw std::runtime_error("missing residue");
                    }
                    v[i]=it->second;
                }
            }
        }
    };

}}}}

namespace desres { namespace msys { namespace atomsel {
    KeywordPtr keyword_fragment( SystemPtr ent ) {
        return KeywordPtr(new Fragment(ent));
    }
}}}

