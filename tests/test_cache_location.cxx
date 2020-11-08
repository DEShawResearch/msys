#include <msys/molfile/dtrplugin.hxx>

using namespace desres::molfile;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        auto loc = StkReader::filename_to_cache_location_v8(argv[i]);
        printf("%s\t\t%s\n", argv[i], loc.data());
    }
    return 0;
}
