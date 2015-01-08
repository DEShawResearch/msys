#include "archive.hxx"
#include "io.hxx"
#include <stdio.h>

using namespace desres::msys;

static void inout(SystemPtr mol) {
    std::stringstream ss;
    SaveArchive(mol, ss);
    SystemPtr newmol = LoadArchive(ss);
    printf("got newmol %p with %lu tables\n", 
            newmol.get(), newmol->tableNames().size());
}
int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        SystemPtr mol = Load(argv[i]);
        inout(mol);
    }

    /* shared param tables */
    SystemPtr mol = System::create();
    ParamTablePtr p = ParamTable::create();
    mol->addTable("a", 1, p);
    mol->addTable("b", 1, p);
    inout(mol);

    return 0;
}

