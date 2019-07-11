#include "io.hxx"
#include "json.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        auto m = ImportJson(argv[i]);
        //std::cout << FormatJson(m, Provenance::fromArgs(argc, argv)) << '\n';
    }
    return 0;
}


