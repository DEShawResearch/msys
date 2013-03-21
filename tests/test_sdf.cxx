#include "sdf.hxx"
#include <stdio.h>
#include <iostream>
#include <fstream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc!=3) {
        fprintf(stderr, "Usage: %s input output.sdf\n", argv[0]);
        return 1;
    }
    SystemPtr mol=Load(argv[1]);
    ExportSdf(mol, argv[2]);
    return 0;
}


