#include "io.hxx"
#include <iostream>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        auto loader = IndexedFileLoader::create(argv[i]);
        std::cout << loader->path() << ": " << loader->size() << " entries.\n";
        for (Id j=loader->size(); j!=0; j--) {
            auto mol = loader->at(j-1);
            std::cout << "entry " << j << " " << mol->name << "\n";
        }
    }
    return 0;
}

