#include "prmtop.hxx"
#include "dms.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        printf("importing %s\n", argv[i]);
        SystemPtr mol = ImportPrmTop(argv[i]);
        printf("exporting %s\n", argv[i]);
        ExportDMS(mol, "prm.dms", Provenance::fromArgs(argc, argv));
    }
    return 0;
}
