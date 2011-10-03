#include "dms.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc>1) {
        SystemPtr sys = ImportDMS(argv[1]);
        sys->renumberGids();
        if (argc>2) {
            ExportDMS(sys, argv[2]);
        }
    }
    return 0;
}

