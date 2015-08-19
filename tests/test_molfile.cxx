#include "molfile/molfile.hxx"
#include <stdio.h>

namespace mf = desres::molfile;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        const char* path = argv[i];
        const molfile_plugin_t* plugin = mf::plugin_for_path(path);
        if (!plugin) {
            fprintf(stderr, "No plugin for '%s'\n", path);
            continue;
        }
        mf::Reader r(plugin, path);
        printf("Got Reader for %s: natoms %ld nframes %ld\n", 
                path, r.natoms(), r.nframes());
        for (;;) {
            auto f = r.next();
            if (!f) break;
            printf("got frame with time %g\n", f->time());
            break;
        }
    }
    return 0;
}
