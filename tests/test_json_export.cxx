#include "io.hxx"
#include <iostream>
#include "json.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        FileFormat fmt;
        SystemPtr m = Load(argv[i], false, &fmt);
        if (!m) { 
            MSYS_FAIL("Unable to guess filetype of " << argv[i]);
        }
        std::cout << FormatJson(m, Provenance::fromArgs(argc, argv)) << '\n';
    }
    return 0;
}


