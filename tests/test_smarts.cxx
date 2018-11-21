#include "smarts.hxx"
#include "io.hxx"
#include "analyze.hxx"
#include "atomsel.hxx"

#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<3) {
        fprintf(stderr, "Usage: %s input.dms smarts1 [smarts2...]\n", argv[0]);
        return 1;
    }
    SystemPtr mol = Load(argv[1]);
    AssignBondOrderAndFormalCharge(mol);
    std::unique_ptr<AnnotatedSystem> annot_mol(new AnnotatedSystem(mol));
    IdList sel = mol->atoms();
    for (int i=2; i<argc; i++) {
        double t=-now();
        MultiIdList matches=SmartsPattern(argv[i]).findMatches(*annot_mol, sel);
        t+=now();
        printf("%s: %lu matches in %.3fms\n", argv[i], matches.size(), t*1000);
    }
    return 0;
}
