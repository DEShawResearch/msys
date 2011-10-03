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

    bool is_water( SystemPtr sys, Id res ) {
        Id O(BadId), H1(BadId), H2(BadId);
        IdList atoms = sys->atomsForResidue(res);
        for (Id i=0; i<atoms.size(); i++) {
            Id id=atoms[i];
            const atom_t& b = sys->atom(id);
            if (b.atomic_number==8) {
                if (bad(O)) O=id;
                else {
                    O=BadId;
                    break;
                }
            } else if (b.atomic_number==1) {
                if      (bad(H1)) H1=id;
                else if (bad(H2)) H2=id;
                else {
                    O=BadId;
                    break;
                }
            }
        }
        if (bad(O) || bad(H1) || bad(H2)) return false;
        return
            sys->bondCountForAtom(O)==2 &&
            sys->bondCountForAtom(H1)==1 &&
            sys->bondCountForAtom(H2)==1 &&
            sys->bond(sys->bondsForAtom(H1)[0]).other(H1)==O &&
            sys->bond(sys->bondsForAtom(H2)[0]).other(H2)==O;
    }

    bool has_water_residue_name( const std::string& resname ) {
        static const char * names[] = {
            "H2O", "HH0", "OHH", "HOH", "OH2", "SOL", "WAT", "TIP",
            "TIP2", "TIP3", "TIP4", "SPC"
        };
        unsigned i,n = sizeof(names)/sizeof(names[0]);
        for (i=0; i<n; i++) {
            if (resname==names[i]) return true;
        }
        return false;
    }

    class Water : public Keyword {
        SystemPtr _sys;
    public:
        Water(SystemPtr sys) : Keyword("water", KEY_INT), _sys(sys) {}
        void iget(const Selection& s, std::vector<Int>& v) const {
            typedef std::map<Id,bool> ResMap;
            ResMap resmap;
            for (Id i=0; i<s.size(); i++) {
                if (!s[i]) continue;
                Id res = _sys->atom(i).residue;
                std::pair<ResMap::iterator,bool> r=resmap.insert(
                        ResMap::value_type(res,false));
                if (r.second) {
                    const std::string& resname = _sys->residue(res).name;
                    /* res was a new value, so do water check */
                    r.first->second=
                        is_water(_sys,res) || has_water_residue_name(resname);
                }
                v[i]=r.first->second;
            }
        }
    };

    enum BackboneType { ATOM_NORMAL, ATOM_PBACK, ATOM_NBACK };
    enum ResidueType { RESIDUE_NOTHING, RESIDUE_PROTEIN, RESIDUE_NUCLEIC };
    typedef std::map<Id, BackboneType> AtypeMap;

    ResidueType analyze_residue( SystemPtr sys, Id res, AtypeMap& amap ) { 

        static const char * protypes[] = { "CA", "C", "O", "N" };
        /* NOTE: VMD 1.7 _tries_ to include O2 in its backbone selection, but
         * because of a bug in the implementation, it doesn't include it.  We
         * match that behavior here. */
        static const char * proterms[] = {
            "OT1", "OT2", "OXT", "O1" /* , "O2" */
        };
        static const char * nuctypes[] = {
            "P", "O1P", "O2P", "OP1", "OP2", "C3*", "C3'", "O3*", "O3'",
            "C4*", "C4'", "C5*", "C5'", "O5*", "O5'"
        };
        static const char * nucterms[] = {
            "H5T", "H3T"
        };
        typedef std::map<std::string,BackboneType> NameMap;
        static NameMap types, terms;
        if (!types.size()) {
            for (unsigned i=0; i<sizeof(protypes)/sizeof(protypes[0]); i++) {
                types[protypes[i]]=ATOM_PBACK;
            }
            for (unsigned i=0; i<sizeof(nuctypes)/sizeof(nuctypes[0]); i++) {
                types[nuctypes[i]]=ATOM_NBACK;
            }
            for (unsigned i=0; i<sizeof(proterms)/sizeof(proterms[0]); i++) {
                terms[proterms[i]]=ATOM_PBACK;
            }
            for (unsigned i=0; i<sizeof(nucterms)/sizeof(nucterms[0]); i++) {
                terms[nucterms[i]]=ATOM_NBACK;
            }
        }

        int npro=0, nnuc=0;
        IdList atoms = sys->atomsForResidue(res);
        for (Id i=0; i<atoms.size(); i++) {
            Id id = atoms[i];
            const atom_t& atm = sys->atom(id);
            /* ignore pseudos in determination of residue type */
            if (!atm.atomic_number) continue;
            const std::string& aname = atm.name;
            /* check for nucleic or protein backbone */
            NameMap::const_iterator iter=types.find(aname);
            BackboneType atype=ATOM_NORMAL;
            if (iter!=types.end()) {
                atype=iter->second;
            } else {
                /* try terminal names */
                iter=terms.find(aname);
                if (iter!=terms.end()) {
                    /* must be bonded to atom of the same type */
                    IdList bonds = sys->bondsForAtom(id);
                    for (Id j=0; j<bonds.size(); j++) {
                        Id o = sys->bond(bonds[j]).other(id);
                        AtypeMap::const_iterator oiter=amap.find(o);
                        if (oiter!=amap.end() && oiter->second==iter->second) {
                            atype=iter->second;
                            break;
                        }
                    }
                }
            }
            amap[id]=atype;
            if (atype==ATOM_PBACK) ++npro;
            if (atype==ATOM_NBACK) ++nnuc;
        }
        ResidueType rtype=RESIDUE_NOTHING;
        if      (npro>=4) rtype=RESIDUE_PROTEIN;
        else if (nnuc>=4) rtype=RESIDUE_NUCLEIC;
        else {
            /* make the atom types consistent with the residue type.  Note
             * that this isn't quite right - it's possible that we could
             * have atoms labeled as PBACK in a NUCLEIC residue, or vice
             * versa, but this is what VMD does so we will, too.  */
            for (AtypeMap::iterator i=amap.begin(); i!=amap.end(); ++i) 
                i->second=ATOM_NORMAL;
        }
        return rtype;
    }
    class Backbone : public Keyword {
        SystemPtr _sys;
    public:
        Backbone(SystemPtr sys)
        : Keyword("backbone", KEY_INT), _sys(sys) {}
        void iget(const Selection& s, std::vector<Int>& v) const {
            typedef std::map<Id,bool> AtmMap;
            AtmMap atmmap;
            for (Id atm=0; atm<s.size(); atm++) {
                if (!s[atm]) continue;
                if (!atmmap.count(atm)) {
                    /* haven't seen this residue before; analyze it */
                    AtypeMap atypes;
                    analyze_residue(_sys, _sys->atom(atm).residue, atypes);
                    for (AtypeMap::const_iterator iter=atypes.begin(); 
                            iter!=atypes.end(); ++iter) 
                        atmmap[iter->first] = 
                            iter->second==ATOM_PBACK || iter->second==ATOM_NBACK;
                }
                v[atm]=atmmap[atm];
            }
        }
    };
    class ResType : public Keyword {
        SystemPtr _sys;
        ResidueType _type;
        public:
        ResType(SystemPtr sys, ResidueType type) 
            : Keyword( type==RESIDUE_PROTEIN ? "protein" :
                    type==RESIDUE_NUCLEIC ? "nucleic" :
                    "nothing", KEY_INT ),
            _sys(sys),
            _type(type) {}
        void iget(const Selection& s, std::vector<Int>& v) const {
            typedef std::map<Id,bool> ResMap;
            ResMap resmap;
            for (Id i=0; i<s.size(); i++) {
                if (!s[i]) continue;
                Id res = _sys->atom(i).residue;
                std::pair<ResMap::iterator,bool> r=resmap.insert(
                        ResMap::value_type(res,false));
                if (r.second) {
                    /* res was a new value, so do residue analysis */
                    AtypeMap amap;
                    ResidueType rtype=analyze_residue(_sys,res,amap);
                    r.first->second = (rtype==_type);
                }
                v[i]=r.first->second;
            }
        }
    };

}}}}

namespace desres { namespace msys { namespace atomsel {
    KeywordPtr keyword_fragment( SystemPtr ent ) {
        return KeywordPtr(new Fragment(ent));
    }
    KeywordPtr keyword_water( SystemPtr ent ) {
        return KeywordPtr(new Water(ent));
    }
    KeywordPtr keyword_backbone( SystemPtr ent ) {
        return KeywordPtr(new Backbone(ent));
    }
    KeywordPtr keyword_protein( SystemPtr ent ) {
        return KeywordPtr(new ResType(ent, RESIDUE_PROTEIN));
    }
    KeywordPtr keyword_nucleic( SystemPtr ent ) {
        return KeywordPtr(new ResType(ent, RESIDUE_NUCLEIC));
    }
}}}

