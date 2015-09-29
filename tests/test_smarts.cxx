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
    //AnnotatedSystem annot_mol(mol);
    std::unique_ptr<AnnotatedSystem> annot_mol(new AnnotatedSystem(mol));
    IdList sel = Atomselect(mol, "not water");
    //IdList sel = mol->atoms();
    for (int i=2; i<argc; i++) {
        printf("%s\n", argv[i]);
        MultiIdList matches=SmartsPattern(argv[i]).findMatches(*annot_mol, sel);
        BOOST_FOREACH(IdList l, matches){
           printf(" --> [ ");
           BOOST_FOREACH(Id i, l){
             printf("%u ",i);
           }
           printf("]\n");
        }
    }
    return 0;
}
