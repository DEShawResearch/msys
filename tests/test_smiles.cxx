#include "smiles.hxx"
#include "sdf.hxx"
#include "atomsel.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        const char* s = argv[i];
        auto mol = FromSmilesString(s);
        if (!mol) {
            printf("FAILED: %s\n", s);
            continue;
        }
        mol->name() = s;
        std::cout << FormatSdf(*mol);
    }
    return 0;
}
