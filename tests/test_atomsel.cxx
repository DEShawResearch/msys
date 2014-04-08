#include "dms.hxx"
#include "atomsel.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<3) {
        fprintf(stderr, "Usage: %s input.dms atomsel [atomsel...]\n", argv[0]);
        return 1;
    }
    SystemPtr sys = ImportDMS(argv[1]);
    for (int i=2; i<argc; i++) {
        double t=-now();
        IdList atoms = Atomselect(sys, argv[i]);
        t+=now();
        printf("%s -> %d atoms (%8.3fms)\n", argv[i], (int)atoms.size(), t*1000);
    }
    return 0;
}

