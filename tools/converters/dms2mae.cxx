#include "mae.hxx"
#include "dms.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc!=3) {
        fprintf(stderr, "Usage: %s input.dms output.mae\n", argv[0]);
        exit(1);
    }
    SystemPtr h = ImportDMS(argv[1]);
    ExportMAE(h, argv[2]);
    return 0;
}
