#include "load.hxx"
#include <iostream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        FileFormat fmt;
        SystemPtr m = Load(argv[i], &fmt);
        std::cout << "Guess file format " << FileFormatName(fmt) << std::endl;
        if (!m) { 
            MSYS_FAIL("Unable to guess filetype of " << argv[i]);
        }
    }
    return 0;
}

