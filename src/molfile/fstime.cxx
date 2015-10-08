#include "molfile/dtrplugin.hxx"
#include <stdio.h>

int main(int argc, char *argv[]) {
    for (int i=1; i<argc; i++) {
        desres::molfile::DtrReader r(argv[i]);
        r.init();
        unsigned long size = r.size();
        double last_time=0;
        if (size>0) {
            r.times(size-1, 1, &last_time);
        }
        printf("%lu %.17g\n", size, last_time);
    }
    return 0;
}
