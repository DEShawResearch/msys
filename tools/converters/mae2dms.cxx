#include "mae.hxx"
#include <dms.hxx>
#include <boost/filesystem.hpp>

#include <string.h>
#include <fstream>

using namespace desres::msys;

static bool ignore_unrecognized = false;

static void usage(FILE *fd) {
    fprintf(fd,
"**mae2dms input.mae output.dms**\n"
"\n"
"  Converts an mae file to a dms file, preserving as much forcefield\n"
"  information as possible.\n"
"\n"
"Options:\n"
"   -h, --help              show this help message and exit\n"
"   --ignore-unrecognized   skip unrecognized ffio_ff subblocks\n"
"\n"
"mae2dms converts mae files to dms files.  Atom order and forcefield\n"
"information are all preserved, though there are likely to be differences\n"
"in the precise values of the forces computed from the file due differences\n"
"in how force terms are represented.\n"
);
}


void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "-h") ||
             !strcmp(pargv[0][j], "--help")) {
            usage(stdout);
            exit(0);
        } else if ( !strcmp(pargv[0][j], "--ignore-unrecognized")) {
            ignore_unrecognized = true; 
        }  else {
            j++;
        }
    }
    pargc[0]=j;
    pargv[0][j]=NULL;
}


int main(int argc, char *argv[]) {

    parse_cmdline(&argc, &argv);
    if (argc!=3) {
        usage(stderr);
        exit(1);
    }

    SystemPtr h;
    try {
        h = ImportMAE(argv[1], ignore_unrecognized);
    }
    catch (std::exception& e) {
        fprintf(stderr, "mae2dms: Failed importing file at '%s'\n", argv[1]);
        fprintf(stderr, "ERROR: %s\n", e.what());
        exit(1);
    }


    try {
        ExportDMS(h, argv[2], Provenance::fromArgs(argc, argv));
    }
    catch (std::exception& e) {
        fprintf(stderr, "mae2dms: Failed exporting dms to '%s'\n", argv[2]);
        fprintf(stderr, "ERROR: %s\n", e.what());
        exit(1);
    }
    return 0;
}

