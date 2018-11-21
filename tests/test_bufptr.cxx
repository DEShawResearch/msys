#include "molfile/dtrplugin.hxx"

using namespace desres::molfile;

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        DtrReader r(argv[i]);
        r.init();
        printf("%s: %ld frames\n", argv[i], r.size());
        void* bufptr = nullptr;
        molfile_timestep_t ts[1];
        ts->coords = new float[r.natoms()];
        for (ssize_t fid=0, n=r.size(); fid<n; fid++) {
            auto keymap = r.frame(fid, nullptr, &bufptr);
            ///printf("fid %ld keymap size %lu\n", fid, keymap.size());
        }
        free(bufptr);
        delete[] ts->coords;
    }
    return 0;
}
