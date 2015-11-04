#include "io.hxx"
#include "atomsel.hxx"
#include "spatial_hash.hxx"
#include <stdio.h>
#include <numeric>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<5) {
        fprintf(stderr, "Usage: %s input.dms target query [rad rad rad...]\n", argv[0]);
        return 1;
    }
    SystemPtr sys = Load(argv[1]);
    std::string target(argv[2]);
    std::string query(argv[3]);

    /* fetch all positions */
    std::vector<float> pos;
    sys->getPositions(std::back_inserter(pos));
    const double* cell = sys->global_cell[0];

    /* evaluate selections */
    IdList tsel = Atomselect(sys, target);
    IdList qsel = Atomselect(sys, query);
    printf("tsel %lu qsel %lu\n", tsel.size(), qsel.size());

    /* extract positions into contiguous arrays, as in zendo. */
    std::vector<float> tpos, qpos;
    for (auto id : tsel) {
        const float* p = &pos[3*id];
        tpos.insert(tpos.end(), p, p+3);
    }
    for (auto id : qsel) {
        const float* p = &pos[3*id];
        qpos.insert(qpos.end(), p, p+3);
    }

    /* need to generate bogus ids */
    IdList tids(tsel.size()), qids(qsel.size());
    std::iota(tids.begin(), tids.end(), 0);
    std::iota(qids.begin(), qids.end(), 0);

    /* result of selection */

    for (int i=4; i<argc; i++) {
        double t;
        double rad = atof(argv[i]);

        {
        t=-now();
        auto S = FindWithin(qids, &qpos[0], tids, &tpos[0], rad, NULL);
        t+=now();
        printf("within   -> %lu atoms (%8.3fms)\n", S.size(), t*1000);
        }

        if (0) {
        t=-now();
        auto S = FindWithin(qids, &qpos[0], tids, &tpos[0], rad, cell);
        t+=now();
        printf("pbwithin -> %lu atoms (%8.3fms)\n", S.size(), t*1000);
        }

    }
    return 0;
}

