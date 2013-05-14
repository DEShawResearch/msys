#include "mae.hxx"
#include "io.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc>1) {
        SystemPtr sys = Load(argv[1]);
        if (argc>2) {
            ExportMAE(sys, argv[2], Provenance::fromArgs(argc,argv));
        }
    }
    return 0;
}

