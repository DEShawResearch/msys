#include "load.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) Load(argv[i]);
    return 0;
}

