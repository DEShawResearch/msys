#include "mae.hxx"
#include "dms.hxx"
#include <cstdio>
#include <string.h>

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
"   -a, --assign-type       Reassign m_mmod_type\n"
"\n"
"dms2mae is the inverse of mae2dms, though with a more restricted\n"
"range of supported forcefield terms.\n"
);
}

static bool assign_type = false;

static
void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "-h") ||
             !strcmp(pargv[0][j], "--help")) {
            usage(stdout);
            exit(0);
        } else if ( !strcmp(pargv[0][j], "--assign-type") ||
                    !strcmp(pargv[0][j], "-a")) {
            assign_type = true;
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
    SystemPtr h = ImportDMS(argv[1]);
    if (assign_type) h->delAtomProp(h->atomPropIndex("m_mmod_type"));
    ExportMAE(h, argv[2], Provenance::fromArgs(argc,argv));
    return 0;
}
