#include <msys/molfile/dtrplugin.hxx>
#include <msys/system.hxx>
#include <fstream>

using namespace desres::msys;
using namespace desres::molfile;

int main() {
    auto stk_path = "tests/files/run.stk";
    StkReader stk(stk_path);
    stk.init();

    for (Id i=0, n=stk.nframesets(); i<n; i++) {
        auto dtr = stk.frameset(i);
        printf("%8d %s\n", i, dtr->path().data());
        dtr->set_path("/mafs/12345/" + dtr->path().substr(20));
    }

    auto cache_path = "updated.cache";
    stk.write_cachefile(cache_path);

     StkReader stk2(stk_path);
     const bool verbose = true;
     const bool v8 = true;
     stk2.read_stk_cache_file(cache_path, verbose, v8);

    for (Id i=0, n=stk2.nframesets(); i<n; i++) {
        auto dtr = stk2.frameset(i);
        printf("%8d %s\n", i, dtr->path().data());
    }

    return 0;
}

