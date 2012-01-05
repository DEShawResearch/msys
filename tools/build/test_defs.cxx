#include "defs.hxx"

using namespace desres::msys::builder;

int main(int argc, char *argv[]) {

    defs_t defs;

    for (int i=1; i<argc; i++) {
        defs.import_charmm_topology(argv[i]);
    }

    printf("have %d types, %d residues\n", (int)defs.types.size(),
            (int)defs.residues.size());
    return 0;
}

