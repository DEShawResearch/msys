#include "k-nearest.hxx"
#include "within_predicate.hxx"
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
}

void KNearestPredicate::eval( Selection& S ) {

    /* evaluate subselection */
    Selection subsel = full_selection(_sys);
    _sub->eval(subsel);

    /* empty subselection -> empty selection */
    if (!subsel.size()) {
        S.clear();
        return;
    }
    const Id subcnt = subsel.count();
    //printf("sub cnt: %u\n", subcnt);

    /* S.count() <= k: nothing to do */
    const Id cnt = S.count();
    //printf("S count: %u\n", cnt);
    if (cnt<=_N) return;

    double rmin=0;
    Selection smin(S);
    exwithin_predicate(_sys, rmin, _sub)->eval(smin);
    Id nmin = smin.count();

    double rmax=2.5;
    Selection smax(0);
    Id nmax=nmin;

    /* increase rmax until Ny is at least _N */
    for (;;) {
        smax = S;
        exwithin_predicate(_sys, rmax, _sub)->eval(smax);
        nmax = smax.count();
        //printf("rmin %f nmin %u rmax %f nmax %u\n", rmin, nmin, rmax, nmax);
        if (nmax >= _N) break;
        rmin = rmax;
        smin = smax;
        nmin = nmax;
        rmax *= 1.5;
    }

    /* Do a couple rounds of bisection search to narrow it down */
    for (int nb=0; nb<6; nb++) {
        Selection sm(S);
        double rm = 0.5*(rmin+rmax);
        exwithin_predicate(_sys, rm, _sub)->eval(sm);
        Id nm = sm.count();
        //printf("rm %f nm %u\n", rm, nm);
        if (nm>=_N) {
            smax = sm;
            rmax = rm;
            nmax = nm;
        } else {
            smin = sm;
            rmin = rm;
            nmin = nm;
        }
    }
    //printf("min: rad %f n %u\n", rmin, nmin);
    //printf("max: rad %f n %u\n", rmax, nmax);

    /* cache the protein positions.
     * FIXME: consider transpose for SIMD */
    std::vector<double> pro;
    pro.reserve(3*subcnt);
    for (Id i=0, n=S.size(); i<n; i++) {
        if (subsel[i]) {
            const atom_t& atm = _sys->atomFAST(i);
            pro.insert(pro.end(), &atm.x, &atm.x+3);
        }
    }

    //printf("cached %lu pro positions\n", pro.size()/3);

    std::vector<std::pair<double,Id> > pts;
    /* for each water in smax but not in smin */
    for (Id i=0, n=S.size(); i<n; i++) {
        if (smin[i] || !smax[i]) continue;
        double r2 = std::numeric_limits<double>::max();
        double x = _sys->atomFAST(i).x;
        double y = _sys->atomFAST(i).y;
        double z = _sys->atomFAST(i).z;
        /* find min dist to protein */
        const double* p = &pro[0];
        for (Id j=0, m=subcnt; j<m; j++) {
            double dx = x-p[0];
            double dy = y-p[1];
            double dz = z-p[2];
            double d2 = dx*dx + dy*dy + dz*dz;
            r2 = std::min(r2, d2);
            p += 3;
        }
        pts.push_back(std::make_pair(r2, i));
    }
    std::partial_sort(pts.begin(), pts.begin()+(_N-nmin), pts.end());
    S = smin;
    for (Id i=0, n=_N-nmin; i<n; i++) {
        S[pts[i].second] = 1;
    }
}

namespace desres { namespace msys { namespace atomsel {

    PredicatePtr k_nearest_predicate(SystemPtr sys, unsigned k, PredicatePtr S) {
        return PredicatePtr(new KNearestPredicate(sys,k,S));
    }

}}}


