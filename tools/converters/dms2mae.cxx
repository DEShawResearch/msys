#include "mae.hxx"
#include "dms.hxx"
#include <cstdio>

using namespace desres::msys;

static void usage(FILE *fd) {
    fprintf(fd,
"**mae2dms input.dms output.mae**\n"
"\n"
"  Converts a DMS file to an MAE file, preserving as much forcefield\n"
"  information as possible.\n"
"\n"
"Options:\n"
"   -h, --help              show this help message and exit\n"
"\n"
"dms2mae is the inverse of mae2dms, though with a more restricted\n"
"range of supported forcefield terms.\n"
);
}

int main(int argc, char *argv[]) {
    if (argc==2 && argv[1][0]=='-' && argv[1][1]=='h') {
        usage(stdout);
        exit(0);
    }
    if (argc!=3) {
        usage(stderr);
        exit(1);
    }
    SystemPtr h = ImportDMS(argv[1]);
    ExportMAE(h, argv[2], Provenance::fromArgs(argc,argv));
    return 0;
}
