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

namespace desres { namespace msys { namespace atomsel {

    Selection full_selection(SystemPtr sys) {
        if (!sys->atomCount()) return Selection(0,IdList());
        Id last = sys->atoms().back();
        return Selection(last+1, sys->atoms());
    }

}}}

namespace {
  using namespace desres::msys::atomsel;
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

INT_KEY(fragid,sys->atom(i).fragid)
INT_KEY(anum,sys->atom(i).atomic_number)
INT_KEY(hydrogen,sys->atom(i).atomic_number==1)
INT_KEY(numbonds,sys->bondCountForAtom(i))
INT_KEY(resid,sys->residue(sys->atom(i).residue).resid)
INT_KEY(residue,sys->atom(i).residue)
INT_KEY(index,i)
INT_KEY(backbone,sys->atom(i).type==AtomProBack || sys->atom(i).type==AtomNucBack)
INT_KEY(water,sys->residue(sys->atom(i).residue).type==ResidueWater)
INT_KEY(protein,sys->residue(sys->atom(i).residue).type==ResidueProtein)
INT_KEY(nucleic,sys->residue(sys->atom(i).residue).type==ResidueNucleic)


DBL_KEY(charge,atom(i).charge)
DBL_KEY(mass,atom(i).mass)
DBL_KEY(x,atom(i).x)
DBL_KEY(y,atom(i).y)
DBL_KEY(z,atom(i).z)
DBL_KEY(vx,atom(i).vx)
DBL_KEY(vy,atom(i).vy)
DBL_KEY(vz,atom(i).vz)

STR_KEY(chain,sys->chain(sys->residue(sys->atom(i).residue).chain).name)
STR_KEY(name,sys->atom(i).name)
STR_KEY(resname,sys->residue(sys->atom(i).residue).name)
STR_KEY(element,desres::msys::AbbreviationForElement(sys->atom(i).atomic_number))

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


