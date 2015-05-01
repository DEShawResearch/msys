#include "smiles.hxx"
#include "io.hxx"
#include "atomsel.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        const char* s = argv[i];
        printf("%s\n", s);
        SystemPtr mol = FromSmilesString(s);
        if (!mol) {
            printf("FAILED\n");
            continue;
        }
        Save(mol, "stdout.sdf", Provenance::fromArgs(argc,argv), 0);
    }
    return 0;
}
