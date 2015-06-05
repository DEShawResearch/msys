#include "io.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        //printf("reading %s\n", argv[i]);
        LoadIteratorPtr iter = LoadIterator::create(argv[i]);
        SystemPtr mol;
        while ((mol=iter->next())) {
            //printf("  molecule '%s': %u atoms, %u bonds\n",
                    //mol->name.c_str(), mol->atomCount(), mol->bondCount());
        }
    }
    return 0;
}

