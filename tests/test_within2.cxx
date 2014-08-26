#include "io.hxx"
#include "atomsel.hxx"
#include "atomsel/within_predicate.hxx"
#include "within.hxx"

#include <stdio.h>
#include <boost/foreach.hpp>


using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<5) {
        fprintf(stderr, "Usage: %s input.dms water protein radius*\n", argv[0]);
        return 1;
    }
    SystemPtr sys = Load(argv[1]);
    IdList water = Atomselect(sys, argv[2]);
    IdList protein = Atomselect(sys, argv[3]);
    std::vector<float> coords(sys->atomCount()*3);
    sys->getPositions(coords.begin());
    const float* pos = &coords[0];

    for (int i=4; i<argc; i++) {
        float radius = atof(argv[i]);

        double t2=-now();
        unsigned ns = SpatialHash<float>(pos, protein)
            .findWithin(radius, pos, water).size();
        t2 += now();

        double t=-now();
        IdList sel = atomsel::FindWithin(water, pos, protein, pos, radius, 0);
        t += now();

        printf("FindWithin  %5.2f -> %6lu atoms (%8.3fms) atoms/ms %8.3f\n", radius, 
                sel.size(), 1000*t, sel.size()/(1000*t));

        printf("SpatialHash %5.2f -> %6u atoms (%8.3fms) atoms/ms %8.3f\n", radius, 
                ns, 1000*t2, ns/(1000*t2));

    }
    return 0;
}

