#include "dtrplugin.hxx"
#include <stdio.h>

static void print_usage() {
    printf("Usage: fstime input.dtr [input.dtr ...]\n");
}

int main(int argc, char *argv[]) {
    if (argc==1 || (argc==2 && (!strcmp(argv[1], "-h") ||
                                !strcmp(argv[1], "--help")))) {
        print_usage();
        return 0;
    }
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
