#include "amber.hxx"
#include "dms.hxx"
#include <sys/stat.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        printf("importing %s\n", argv[i]);
        SystemPtr mol = ImportPrmTop(argv[i]);
        std::string crd(argv[i]);
        crd.resize(crd.size()-6);
        crd += "crd";
        struct stat buf[1];
        if (!stat(crd.c_str(), buf)) {
            printf("reading crd %s\n", crd.c_str());
            ImportCrdCoordinates(mol, crd);
        }
        printf("exporting %s\n", argv[i]);
        ExportDMS(mol, "prm.dms", Provenance::fromArgs(argc, argv));
    }
    return 0;
}
