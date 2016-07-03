#include "smiles.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        auto mol = FromSmilesString(argv[i]);
        printf("%s: %u atoms %u bonds\n", argv[i], mol->atomCount(), mol->bondCount());
    }
    return 0;
}

