#include "dms.hxx"
using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        std::cout << argv[i] << " " << HashDMS(argv[i]) << "\n";
    }
    return 0;
}
