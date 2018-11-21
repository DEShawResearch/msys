#include "io.hxx"
#include <iostream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<2) exit(1);
    auto loader = IndexedFileLoader::create(argv[1]);
    std::cout << loader->path() << ": " << loader->size() << " entries.\n";
    for (int i=2; i<argc; i++) {
        size_t index = atoi(argv[i]);
        auto mol = loader->at(index);
        std::cout << "entry " << index << " " << mol->name << " " << mol->atomCount() << " atoms\n";
    }
    return 0;
}

