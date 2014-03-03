#include <stdio.h>
#include "cealign.hxx"
#include "io.hxx"
#include "atomsel.hxx"

using namespace desres;
using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<3) {
        fprintf(stderr, "Usage: %s reference.dms system.dms [ output.dms ]\n",
                argv[0]);
        return 1;
    }
    printf("Reading %s\n", argv[1]);
    SystemPtr m1 = Load(argv[1]); 
    printf("Reading %s\n", argv[2]);
    SystemPtr m2 = Load(argv[2]);
    IdList ca1 = Atomselect(m1, "protein and name CA");
    IdList ca2 = Atomselect(m2, "protein and name CA");
    if (ca1.size()<8 || ca2.size()<8) {
        fprintf(stderr, "Systems contain too few CA atoms\n");
        return 1;
    }
    printf("%s: %lu atoms\n", argv[1], ca1.size());
    printf("%s: %lu atoms\n", argv[2], ca2.size());
    std::vector<double> p1, p2;
    m1->getPositions(std::back_inserter(p1));
    m2->getPositions(std::back_inserter(p2));

#if 0
    CEAlign<double> engine = CEAlign<double>::WithDefaults();
    double R[9], T1[3], T2[3];
    double rmsd = engine.compute(ca1.size(), &ca1[0], &p1[0],
                                 ca2.size(), &ca2[0], &p2[0],
                                 R, T1, T2);
    printf("got rmsd %g\n", rmsd);
    if (argc>3) {
        pfx::apply_shift(m2->atomCount(), &p2[0], -T2[0], -T2[1], -T2[2]);
        pfx::apply_rotation(m2->atomCount(), &p2[0], R);
        pfx::apply_shift(m2->atomCount(), &p2[0], T1[0], T1[1], T1[2]);
        m2->setPositions(p2.begin());
        printf("Saving to %s\n", argv[3]);
        Save(m2, argv[3], Provenance::fromArgs(argc,argv), 0);
    }
#endif
    return 0;
}

