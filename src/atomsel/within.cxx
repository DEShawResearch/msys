/* @COPYRIGHT@ */

#include "token.hxx"
#include "../spatial_hash.hxx"

using namespace desres::msys;
using namespace desres::msys::atomsel;

namespace {
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

