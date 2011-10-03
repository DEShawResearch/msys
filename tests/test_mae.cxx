#include "mae.hxx"
#include <fstream>

int main(int argc, char *argv[]) {
    using namespace desres::msys;
    for (int i=1; i<argc; i++) {
        ImportMAE(argv[i]);
    }
    return 0;
}

