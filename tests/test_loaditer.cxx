#include "io.hxx"
#include "clone.hxx"
#include "sdf.hxx"
#include "analyze.hxx"
#include "annotated_system.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        //printf("reading %s\n", argv[i]);
        LoadIteratorPtr iter = LoadIterator::create(argv[i]);
        SystemPtr mol;
        int n=0;
        while ((mol=iter->next())) {
            ++n;
            auto orig = Clone(mol, mol->atoms());
            try {
                AddHydrogens(mol, mol->atoms());
                AnnotatedSystem::create(mol);
            } catch (std::exception& e) {
                fprintf(stderr, "Failed system %d: %s\n%s\n", 
                        n, mol->name.c_str(), e.what());
                ExportSdf(orig, std::cout);
            }
            //printf("  molecule '%s': %u atoms, %u bonds\n",
                    //mol->name.c_str(), mol->atomCount(), mol->bondCount());
        }
    }
    return 0;
}

