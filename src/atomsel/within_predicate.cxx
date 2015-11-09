/* @COPYRIGHT@ */

#include "within_predicate.hxx"
#include "msys_keyword.hxx"
#include "../spatial_hash.hxx"

using namespace desres::msys;
using namespace desres::msys::atomsel;

namespace {
  class WithinPredicate : public Predicate {
    SystemPtr sys;
    const float* pos;
    const double* cell;
    const float rad;
    PredicatePtr sub;
    const bool exclude;
    const bool periodic;

  public:
    WithinPredicate( SystemPtr e, const float* pos, const double* cell, float r, bool excl, bool per, PredicatePtr s )
      : sys(e), pos(pos), cell(cell), rad(r), sub(s), exclude(excl), periodic(per) {}

    void eval( Selection& s );
  };
  class WithinBondsPredicate : public Predicate {
    SystemPtr sys;
    const int N;
    PredicatePtr sub;

  public:
    WithinBondsPredicate( SystemPtr e, int n, PredicatePtr s )
      : sys(e), N(n), sub(s) {}

    void eval( Selection& s );
  };
  class KNearestPredicate : public Predicate {
    SystemPtr _sys;
    const float* pos;
    const double* cell;
    const unsigned _N;
    const bool periodic;
    PredicatePtr _sub;

  public:
    KNearestPredicate(SystemPtr sys, const float* pos, const double* cell, unsigned k, bool per, PredicatePtr sub)
    : _sys(sys), pos(pos), cell(cell), _N(k), periodic(per), _sub(sub) {}

    void eval(Selection& s);
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
    if (!pos) {
        sys->getPositions(std::back_inserter(coords));
        pos = &coords[0];
    }
    if (periodic && !cell) {
        cell = sys->global_cell[0];
    }
    IdList subsel_ids = subsel.ids();
    IdList S_ids = S.ids();

    IdList ids = FindWithin(S_ids, pos, subsel_ids, pos, rad, cell);
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
    if (!pos) {
        _sys->getPositions(std::back_inserter(coords));
        pos = &coords[0];
    }
    if (periodic && !cell) {
        cell = _sys->global_cell[0];
    }
    IdList subsel_ids = subsel.ids();
    IdList S_ids = S.ids();

    IdList ids = FindNearest(S_ids, pos, subsel_ids, pos, _N, cell);
    S.clear();
    for (Id i=0, n=ids.size(); i<n; i++) S[ids[i]] = 1;
}


namespace desres { namespace msys { namespace atomsel {

  PredicatePtr within_predicate( SystemPtr sys, const float* pos, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,pos,nullptr,r,false,false,s));
  }

  PredicatePtr exwithin_predicate( SystemPtr sys, const float* pos, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,pos,nullptr,r,true,false,s));
  }

  PredicatePtr pbwithin_predicate( SystemPtr sys, const float* pos, const double* cell, double r, PredicatePtr s ) {
    return PredicatePtr(new WithinPredicate(sys,pos,cell,r,false,true,s));
  }

  PredicatePtr withinbonds_predicate( SystemPtr sys, int n, PredicatePtr s ) {
    return PredicatePtr(new WithinBondsPredicate(sys,n,s));
  }

    PredicatePtr k_nearest_predicate(SystemPtr sys, const float* pos, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,pos,nullptr,k,false,S));
    }
    PredicatePtr k_pbnearest_predicate(SystemPtr sys, const float* pos, const double* cell, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,pos,cell,k,true,S));
    }

}}}


