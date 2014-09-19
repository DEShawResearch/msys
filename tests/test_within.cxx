#include "dms.hxx"
#include "atomsel.hxx"
#include "spatial_hash.hxx"
#include <stdio.h>
#include <assert.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<5) {
        fprintf(stderr, "Usage: %s input.dms prosel watsel [rad rad rad...]\n", argv[0]);
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

    IdList psel = Atomselect(sys, prosel);
    IdList wsel = Atomselect(sys, watsel);
    printf("psel %lu wsel %lu\n", psel.size(), wsel.size());

    for (int i=4; i<argc; i++) {
        std::string sel("(");
        sel += watsel + ") and within " + argv[i] + " of (" + prosel + ")";

        double t=-now();
        IdList atoms = Atomselect(sys, sel);
        t+=now();
        printf("%s -> %d atoms (%8.3fms)\n", sel.c_str(), (int)atoms.size(), t*1000);
        double rad = atof(argv[i]);
        t=-now();
        IdList S = FindWithin(wsel, &pos[0], psel, &pos[0], rad, NULL);
        t+=now();
        printf("FindWithin -> %lu atoms (%8.3fms)\n", S.size(), t*1000);
        assert(S==atoms);
    }
    for (int i=4; i<argc; i++) {
        std::string sel("(");
        sel += watsel + ") and pbwithin " + argv[i] + " of (" + prosel + ")";

        double t=-now();
        IdList atoms = Atomselect(sys, sel);
        t+=now();
        printf("%s -> %d atoms (%8.3fms)\n", sel.c_str(), (int)atoms.size(), t*1000);
        double rad = atof(argv[i]);
        t=-now();
        IdList S = FindWithin(wsel, &pos[0], psel, &pos[0], rad, cell);
        t+=now();
        printf("FindWithin -> %lu atoms (%8.3fms)\n", S.size(), t*1000);
        assert(S==atoms);
    }
    return 0;
}

