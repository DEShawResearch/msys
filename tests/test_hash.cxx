#include "io.hxx"
#include "hash.hxx"
#include <iostream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        auto mol = Load(argv[i]);
        std::cout << argv[i] << ' ' << HashSystem(mol) << '\n';
    }
    return 0;
}

