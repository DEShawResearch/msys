#include <dms/dms.hxx>
#include <stdio.h>

using namespace desres::msys;

static void usage(FILE *fd) {
    fprintf(fd,
"Usage: dms-version input.dms\n"
"\n"
"  Writes the DMS version number to stdout.\n"
"\n"
"Options:\n"
"   -h, --help       show this help message and exit\n"
);
}

int main(int argc, char *argv[]) {
    if (argc!=2) {
        usage(stderr);
        exit(1);
    }
    if (argv[1][0]=='-' && argv[1][1]=='h') {
        usage(stdout);
        exit(0);
    }

    const char* path = argv[1];
    Sqlite dms = Sqlite::read(path);
    Reader r = dms.fetch("dms_version");
    if (!r.size()) {
        printf("%s unknown\n", path);
    } else {
        int MAJOR = r.column("major");
        int MINOR = r.column("minor");
        if (MAJOR<0 || MINOR<0) {
            MSYS_FAIL("dms_version table is malformatted in " << path);
        }
        int major = r.get_int(MAJOR);
        int minor = r.get_int(MINOR);
        printf("%s %d.%d\n", path, major, minor);
    }
    return 0;
}
