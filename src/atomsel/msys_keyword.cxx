/* @COPYRIGHT@ */

#include "msys_keyword.hxx"
#include "../elements.hxx"

using desres::msys::SystemPtr;
using desres::msys::Id;
using desres::msys::ResidueWater;
using desres::msys::ResidueProtein;
using desres::msys::ResidueNucleic;
using desres::msys::AtomProBack;
using desres::msys::AtomNucBack;
using desres::msys::AtomProSide;

namespace desres { namespace msys { namespace atomsel {

    Selection full_selection(SystemPtr sys) {
        Selection s(sys->maxAtomId());
        s.fill();
        if (sys->maxAtomId() != sys->atomCount()) {
            for (Id i=0, n=sys->maxAtomId(); i<n; i++) {
                if (!sys->hasAtom(i)) s[i]=0;
            }
        }
        return s;
    }

}}}

namespace {
  using namespace desres::msys::atomsel;

  struct AllKeyword : Keyword {
      AllKeyword() : Keyword("all", KEY_INT) {}
      void iget(const Selection& s, std::vector<Int>& v) const {
        Id i,n=s.size();
        for (i=0; i<n; i++) if (s[i]) v[i]=1;
      }
  };
  struct NoneKeyword : Keyword {
      NoneKeyword() : Keyword("none", KEY_INT) {}
      void iget(const Selection& s, std::vector<Int>& v) const {
        Id i,n=s.size();
        for (i=0; i<n; i++) if (s[i]) v[i]=0;
      }
  };

  struct MsysKeyword : Keyword {
    MsysKeyword( const std::string& n, KeywordType t, SystemPtr _sys)
      : Keyword(n,t), sys(_sys) {}
    SystemPtr sys;
  };

  struct AtomPropKeyword : Keyword {
    AtomPropKeyword( const std::string& n, KeywordType t, SystemPtr _sys)
    : Keyword(n,t), sys(_sys), col(_sys->atomPropIndex(n)) {}

    SystemPtr sys;
    Id col;

    void iget(const Selection& s, std::vector<Int>& v) const {
        Id i,n=s.size();
        for (i=0; i<n; i++) if (s[i]) v[i]=sys->atomPropValue(i,col).asInt();
    }
    void dget(const Selection& s, std::vector<Dbl>& v) const {
        Id i,n=s.size();
        for (i=0; i<n; i++) if (s[i]) v[i]=sys->atomPropValue(i,col).asFloat();
    }
    void sget(const Selection& s, std::vector<Str>& v) const {
        Id i,n=s.size();
        for (i=0; i<n; i++) if (s[i]) v[i]=sys->atomPropValue(i,col).asString();
    }
  };
}

#define INT_KEY(attr,path) \
  namespace { \
    struct Keyword_##attr : MsysKeyword { \
      Keyword_##attr(SystemPtr _ent) : MsysKeyword(#attr,KEY_INT,_ent) {} \
      void iget(const Selection& s, std::vector<Int>& v) const { \
        Id i,n=s.size(); \
        for (i=0; i<n; i++) if (s[i]) v[i]=path; \
      } \
    }; \
  } \
desres::msys::atomsel::KeywordPtr \
desres::msys::atomsel::keyword_##attr( SystemPtr ent ) { \
  return KeywordPtr(new Keyword_##attr(ent)); \
}

#define DBL_KEY(attr,path) \
  namespace { \
    struct Keyword_##attr : MsysKeyword { \
      Keyword_##attr(SystemPtr _ent) : MsysKeyword(#attr,KEY_DBL,_ent) {} \
      void dget(const Selection& s, std::vector<Dbl>& v) const { \
        Id i,n=s.size(); \
        for (i=0; i<n; i++) if (s[i]) v[i]=sys->path; \
      } \
    }; \
  } \
desres::msys::atomsel::KeywordPtr \
desres::msys::atomsel::keyword_##attr( SystemPtr ent ) { \
  return KeywordPtr(new Keyword_##attr(ent)); \
}

#define STR_KEY(attr,path) \
  namespace { \
    struct Keyword_##attr : MsysKeyword { \
      Keyword_##attr(SystemPtr _ent) : MsysKeyword(#attr,KEY_STR,_ent) {} \
      void sget(const Selection& s, std::vector<Str>& v) const { \
        Id i,n=s.size(); \
        for (i=0; i<n; i++) if (s[i]) v[i]=path; \
      } \
    }; \
  } \
desres::msys::atomsel::KeywordPtr \
desres::msys::atomsel::keyword_##attr( SystemPtr ent ) { \
  return KeywordPtr(new Keyword_##attr(ent)); \
}

INT_KEY(fragid,sys->atomFAST(i).fragid)
INT_KEY(anum,sys->atomFAST(i).atomic_number)
INT_KEY(hydrogen,sys->atomFAST(i).atomic_number==1)
INT_KEY(numbonds,sys->bondCountForAtom(i))
INT_KEY(resid,sys->residueFAST(sys->atomFAST(i).residue).resid)
INT_KEY(residue,sys->atomFAST(i).residue)
INT_KEY(index,i)
INT_KEY(backbone,sys->atomFAST(i).type==AtomProBack || sys->atomFAST(i).type==AtomNucBack)
INT_KEY(sidechain,sys->atomFAST(i).type==AtomProSide)
INT_KEY(water,sys->residueFAST(sys->atomFAST(i).residue).type==ResidueWater)
INT_KEY(protein,sys->residueFAST(sys->atomFAST(i).residue).type==ResidueProtein)
INT_KEY(nucleic,sys->residueFAST(sys->atomFAST(i).residue).type==ResidueNucleic)

/* add 1 to ct to match VMD :(  */
INT_KEY(ctnumber,1+sys->chainFAST(sys->residueFAST(sys->atomFAST(i).residue).chain).ct)


DBL_KEY(charge,atomFAST(i).charge)
DBL_KEY(mass,atomFAST(i).mass)
DBL_KEY(x,atomFAST(i).x)
DBL_KEY(y,atomFAST(i).y)
DBL_KEY(z,atomFAST(i).z)
DBL_KEY(vx,atomFAST(i).vx)
DBL_KEY(vy,atomFAST(i).vy)
DBL_KEY(vz,atomFAST(i).vz)

STR_KEY(chain,sys->chainFAST(sys->residueFAST(sys->atomFAST(i).residue).chain).name)
STR_KEY(segid,sys->chainFAST(sys->residueFAST(sys->atomFAST(i).residue).chain).segid)
STR_KEY(name,sys->atomFAST(i).name)
STR_KEY(resname,sys->residueFAST(sys->atomFAST(i).residue).name)
STR_KEY(insertion,sys->residueFAST(sys->atomFAST(i).residue).insertion)
STR_KEY(element,desres::msys::AbbreviationForElement(sys->atomFAST(i).atomic_number))

KeywordPtr desres::msys::atomsel::keyword_atomprop( 
        SystemPtr ent, String const& prop ) {

    Id col = ent->atomPropIndex(prop);
    if (bad(col)) return KeywordPtr();
    ValueType vt = ent->atomPropType(col);
    KeywordType t = vt==IntType ? KEY_INT :
                    vt==FloatType ? KEY_DBL :
                                      KEY_STR ;

    return KeywordPtr(new AtomPropKeyword(prop, t, ent));
}

KeywordPtr desres::msys::atomsel::keyword_all() {
    return KeywordPtr(new AllKeyword);
}
KeywordPtr desres::msys::atomsel::keyword_none() {
    return KeywordPtr(new NoneKeyword);
}
