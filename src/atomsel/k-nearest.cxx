#include "k-nearest.hxx"
#include "msys_keyword.hxx"

#include <algorithm>
#include <vector>
#include <limits>

using namespace desres::msys::atomsel;
using desres::msys::SystemPtr;
using desres::msys::Id;
using desres::msys::atom_t;

namespace {
    class KNearestPredicate : public Predicate {
        SystemPtr _sys;
        const unsigned _N;
        PredicatePtr _sub;

    public:
        KNearestPredicate(SystemPtr sys, unsigned k, PredicatePtr sub)
            : _sys(sys), _N(k), _sub(sub) {}

        void eval(Selection& s);
        void dump( std::ostream& str ) const {
            str << "nearest " << _N << " to [";
            _sub->dump(str);
            str << "]";
        }
    };

    struct PointDistance {
        float d;
        int i;
        PointDistance() {}
        PointDistance(float _d, int _i) : d(_d), i(_i) {}
        bool operator<(const PointDistance& o) const { return d<o.d; }
    };
    struct Position {
        float x,y,z;
        Position() {}
        Position(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
    };
}

/* 
 * Go through each point in S but not in subsel, and find its minimum 
 * distance to any point in subsel.   At the end we sort and take the 
 * first k values.
 *
 * This is a naive implementation with poor scaling.
 */
void KNearestPredicate::eval( Selection& S ) {

    /* evaluate subselection */
    Selection subsel = full_selection(_sys);
    _sub->eval(subsel);

    /* empty subselection -> empty selection */
    if (!subsel.size()) {
        S.clear();
        return;
    }

    /* cache the required positions */
    std::vector<Position> positions(S.size());
    for (Id i=0; i<S.size(); i++) {
        if (S[i] || subsel[i]) {
            const atom_t& atm = _sys->atom(i);
            positions[i].x = atm.x;
            positions[i].y = atm.y;
            positions[i].z = atm.z;
        }
    }

    std::vector<PointDistance> distances;
    for (Id i=0; i<S.size(); i++) {
        if (subsel[i] || !S[i]) continue;
        /* find minimum distance of this point to the subselection */
        const float x = positions[i].x;
        const float y = positions[i].y;
        const float z = positions[i].z;
        float d2=std::numeric_limits<float>::max();
        for (Id j=0; j<subsel.size(); j++) {
            if (!subsel[j]) continue;
            float dx=x-positions[j].x;
            float dy=y-positions[j].y;
            float dz=z-positions[j].z;
            float d2_j = dx*dx + dy*dy + dz*dz;
            if (d2_j<d2) d2=d2_j;
        }
        distances.push_back(PointDistance(d2,i));
    }
    std::sort(distances.begin(), distances.end());

    S.clear();
    unsigned n=_N;
    if (n>distances.size()) n=distances.size();
    for (unsigned i=0; i<n; i++) {
        S[distances[i].i] = 1;
    }
}

namespace desres { namespace msys { namespace atomsel {

    PredicatePtr k_nearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,k,S));
    }

}}}


