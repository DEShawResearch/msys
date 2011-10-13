#include "json.hxx"
#include "dms.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        SystemPtr sys = ImportDMS(argv[i]);
        ExportJSON(sys, std::cout);
        std::cout << std::endl;
    }
    return 0;
}
