#include "load.hxx"
#include "mol2.hxx"
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<2) {
        fprintf(stderr, "Usage: %s input.dms output.mol2\n", argv[0]);
        return 1;
    }
    SystemPtr sys = Load(argv[1]);
    if (argc==3) {
        ExportMol2(sys, argv[2], Provenance::fromArgs(argc,argv));
    }
    return 0;
}

