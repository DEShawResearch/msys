#include "sdf.hxx"
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<2) {
        fprintf(stderr, "Usage: %s input.dms output.mol2\n", argv[0]);
        return 1;
    }
    std::vector<SystemPtr> sdfs= ImportSdfMany(argv[1]);
    return 0;
}


