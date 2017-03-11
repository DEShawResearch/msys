#include "spatial_hash.hxx"
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <cassert>

using namespace desres::msys;

int main() {

    srand48(1973);
    const int N = 10000;
    std::vector<float>  fpos(3*N);
    std::vector<double> dpos(3*N);
    for (int i=0; i<3*N; i++) {
        double d = drand48();
        fpos[i] = d;
        dpos[i] = d;
    }

    IdList A(N/2), B(N/2);
    std::iota(A.begin(), A.end(), 0);
    std::iota(B.begin(), B.end(), N/2);
    double cell[9] = {1,0,0, 0,1,0, 0,0,1};

    SpatialHashT<float> sf(fpos.data(), A.size(), A.data(), cell);
    SpatialHashT<double> sd(dpos.data(), A.size(), A.data(), cell);

    double tf=9999, td=9999;
    for (int i=0; i<10; i++) {
        double t0=now()*1000;
        auto f = sf.findWithin(0.125, fpos.data(), B.size(), B.data());
        double t1=now()*1000;
        auto d = sd.findWithin(0.125, dpos.data(), B.size(), B.data());
        double t2=now()*1000;
        assert(f.size()>0);
        assert(f == d);
        tf=std::min(tf, t1-t0);
        td=std::min(td, t2-t1);
    }
    printf("findWithin: single %.3fms double %.3fms ratio %.3f\n",tf,td,td/tf);

    tf=td=9999;
    for (int i=0; i<10; i++) {
        double t0=now()*1000;
        auto f = sf.findContacts(0.125, fpos.data(), B.size(), B.data());
        double t1=now()*1000;
        auto d = sd.findContacts(0.125, dpos.data(), B.size(), B.data());
        double t2=now()*1000;
        assert(f.size()>0);
        assert(f.size()==d.size());
        for (Id i=0, n=f.size(); i<n; i++) {
            assert(f[i].i == d[i].i);
            assert(f[i].j == d[i].j);
        }
        tf=std::min(tf, t1-t0);
        td=std::min(td, t2-t1);
    }
    printf("findContacts: single %.3fms double %.3fms ratio %.3f\n",tf,td,td/tf);

    return 0;
}

