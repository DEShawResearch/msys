/* @COPYRIGHT@ */

#include "within_predicate.hxx"
#include "msys_keyword.hxx"
#include "../spatial_hash.hxx"

using namespace desres::msys;
using namespace desres::msys::atomsel;

namespace {
  class WithinPredicate : public Predicate {
    SystemPtr sys;
    const float rad;
    PredicatePtr sub;
    const bool exclude;
    const bool periodic;

  public:
    WithinPredicate( SystemPtr e, float r, bool excl, bool per, PredicatePtr s )
      : sys(e), rad(r), sub(s), exclude(excl), periodic(per) {}

    void eval( Selection& s );
    void dump( std::ostream& str ) const {
      if (periodic) str << "pb";
      if (exclude) str << "ex";
      str << "within " << rad << " of [";
      sub->dump(str);
      str << "]";
    }
  };
  class WithinBondsPredicate : public Predicate {
    SystemPtr sys;
    const int N;
    PredicatePtr sub;

  public:
    WithinBondsPredicate( SystemPtr e, int n, PredicatePtr s )
      : sys(e), N(n), sub(s) {}

    void eval( Selection& s );
    void dump( std::ostream& str ) const {
      str << "withinbonds " << N << " of [";
      sub->dump(str);
      str << "]";
    }
  };
  class KNearestPredicate : public Predicate {
    SystemPtr _sys;
    const unsigned _N;
    const bool periodic;
    PredicatePtr _sub;

  public:
    KNearestPredicate(SystemPtr sys, unsigned k, bool per, PredicatePtr sub)
    : _sys(sys), _N(k), periodic(per), _sub(sub) {}

    void eval(Selection& s);
    void dump( std::ostream& str ) const {
        str << (periodic ? "pb" : "")
            << "nearest " << _N << " to [";
        _sub->dump(str);
        str << "]";
    }
  };
}

void WithinPredicate::eval( Selection& S ) {
    Selection subsel = full_selection(sys);
    sub->eval(subsel);
    if (exclude) S.subtract(subsel);

    if (rad<=0) {
        S.intersect(subsel);
        return;
    }

    std::vector<float> coords;
    sys->getPositions(std::back_inserter(coords));
    const float* pos = &coords[0];
    const double* cell = periodic ? sys->global_cell[0] : NULL;
    IdList subsel_ids = subsel.ids();
    IdList S_ids = S.ids();

    IdList ids = SpatialHash(pos, subsel_ids.size(), &subsel_ids[0])
        .findWithin(rad, pos, S_ids.size(), &S_ids[0], cell);

    S.clear();
    for (Id i=0, n=ids.size(); i<n; i++) S[ids[i]] = 1;

}

void WithinBondsPredicate::eval( Selection& S ) {
  Selection subsel = full_selection(sys);
  sub->eval(subsel);
  IdSet atms;
  /* hash the atoms in the subsel */
  for (Id i=0; i<subsel.size(); i++) {
    if (subsel[i]) atms.insert(i);
  }
  /* expand subsel by N bonds */
  for (int i=0; i<N; i++) {
    IdList tmp;
    for (IdSet::const_iterator iter=atms.begin(); iter!=atms.end(); ++iter) {
        IdList bonds = sys->bondsForAtom(*iter);
        for (Id j=0; j<bonds.size(); j++) {
            tmp.push_back(sys->bond(bonds[j]).other(*iter));
        }
    }
    atms.insert(tmp.begin(), tmp.end());
  }
  for (Id i=0; i<subsel.size(); i++) {
    subsel[i] = atms.count(i);
  }
  S.intersect(subsel);
}

void KNearestPredicate::eval( Selection& S ) {

    Selection subsel = full_selection(_sys);
    _sub->eval(subsel);
    S.subtract(subsel);

    std::vector<float> coords;
    _sys->getPositions(std::back_inserter(coords));
    const float* pos = &coords[0];
    const double* cell = periodic ? _sys->global_cell[0] : NULL;
    IdList subsel_ids = subsel.ids();
    IdList S_ids = S.ids();

    IdList ids = SpatialHash(pos, subsel_ids.size(), &subsel_ids[0])
        .findNearest(_N, pos, S_ids.size(), &S_ids[0], cell);

    S.clear();
    for (Id i=0, n=ids.size(); i<n; i++) S[ids[i]] = 1;
}


namespace desres { namespace msys { namespace atomsel {

  PredicatePtr within_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,false,false,s));
  }

  PredicatePtr exwithin_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,true,false,s));
  }

  PredicatePtr pbwithin_predicate( SystemPtr sys, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,r,false,true,s));
  }

  PredicatePtr withinbonds_predicate( SystemPtr sys, int n, PredicatePtr s ) {
    return PredicatePtr(new WithinBondsPredicate(sys,n,s));
  }

    PredicatePtr k_nearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,k,false,S));
    }
    PredicatePtr k_pbnearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,k,true,S));
    }

}}}


