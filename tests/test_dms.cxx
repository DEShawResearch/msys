#include "dms.hxx"
#include "io.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc>1) {
        double t=-now();
        SystemPtr sys = Load(argv[1]);
        t+=now();
        printf("Loaded %s in %8.3fms\n", argv[1], t*1000);
        if (argc>2) {
            ExportDMS(sys, argv[2], Provenance::fromArgs(argc,argv));
        }
    }
    return 0;
}

