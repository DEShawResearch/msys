#include "mae.hxx"
#include <fstream>

int main(int argc, char *argv[]) {
    using namespace desres::msys;
    if (argc==1) {
        SystemPtr mol = ImportMAEFromStream(std::cin);
    } else for (int i=1; i<argc; i++) {
        ImportMAE(argv[i], true);
    }
    return 0;
}

