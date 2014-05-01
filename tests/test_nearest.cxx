#include "dms.hxx"
#include "atomsel.hxx"
#include "atomsel/within_predicate.hxx"
#include <stdio.h>
#include <assert.h>

using namespace desres::msys;
using desres::msys::atomsel::Selection;

int main(int argc, char *argv[]) {
    if (argc<5) {
        fprintf(stderr, "Usage: %s input.dms prosel watsel [k k k...]\n", argv[0]);
        return 1;
    }
    SystemPtr sys = ImportDMS(argv[1]);
    std::string prosel(argv[2]);
    std::string watsel(argv[3]);

    /* fetch all positions */
    std::vector<float> pos;
    sys->getPositions(std::back_inserter(pos));
    const double* cell = sys->global_cell[0];

    /* evaluate selections */
    Selection psel(sys->atomCount());
    Selection wsel(sys->atomCount());
    IdList ids;

    ids = Atomselect(sys, prosel);
    for (Id i=0, n=ids.size(); i<n; i++) psel[ids[i]]=1;
    ids = Atomselect(sys, watsel);
    for (Id i=0, n=ids.size(); i<n; i++) wsel[ids[i]]=1;
    printf("psel %u wsel %u\n", psel.count(), wsel.count());

    for (int i=4; i<argc; i++) {
        std::string sel("(");
        sel += watsel + ") and nearest " + argv[i] + " to (" + prosel + ")";

        double t=-now();
        IdList atoms = Atomselect(sys, sel);
        t+=now();
        printf("%s -> %d atoms (%8.3fms)\n", sel.c_str(), (int)atoms.size(), t*1000);
        Id k= atoi(argv[i]);
        t=-now();
        Selection S = atomsel::FindNearest(wsel, &pos[0], psel, &pos[0],
                                          k, NULL);
        t+=now();
        printf("FindNearest-> %u atoms (%8.3fms)\n", S.count(), t*1000);
        assert(S.ids()==atoms);
    }
    for (int i=4; i<argc; i++) {
        std::string sel("(");
        sel += watsel + ") and pbnearest " + argv[i] + " to (" + prosel + ")";

        double t=-now();
        IdList atoms = Atomselect(sys, sel);
        t+=now();
        printf("%s -> %d atoms (%8.3fms)\n", sel.c_str(), (int)atoms.size(), t*1000);
        Id k = atoi(argv[i]);
        t=-now();
        Selection S = atomsel::FindNearest(wsel, &pos[0], psel, &pos[0],
                                          k, cell);
        t+=now();
        printf("FindNearest -> %u atoms (%8.3fms)\n", S.count(), t*1000);
        assert(S.ids()==atoms);
    }
    return 0;
}

