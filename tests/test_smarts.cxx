#include "smarts.hxx"
#include "load.hxx"
#include "analyze.hxx"
#include "atomsel.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<3) {
        fprintf(stderr, "Usage: %s input.dms smarts1 [smarts2...]\n", argv[0]);
        return 1;
    }
    SystemPtr mol = Load(argv[1]);
    AssignBondOrderAndFormalCharge(mol);
    AnnotatedSystemPtr annot_mol = AnnotatedSystem::create(mol);
    IdList sel = Atomselect(mol, "not water");
    for (int i=2; i<argc; i++) {
        printf("%s\n", argv[i]);
        SmartsPattern(argv[i]).findMatches(annot_mol, sel);
    }
    return 0;
}
