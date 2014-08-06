#include "io.hxx"
#include "analyze.hxx"
#include "mol2.hxx"
#include <cstdio>
#include <string.h>


using namespace desres::msys;

static bool assign_bondorder = true;
static bool assign_sybyltypes = true;

static void usage(FILE *fd) {
    fprintf(fd,
"**msys2mol2 input_file output.mol2**\n"
"\n"
"  Creates a mol2 file from the given input file.\n"
"\n"
"Options:\n"
"   -h, --help       show this help message and exit\n"
"   --preserve-bondorder\n"
"       By default, msys2mol2 will recompute bond orders using the\n"
"       AssignBondOrderAndFormalCharge function.  This option causes\n"
"       msys2mol2 to use whatever bond order information is already\n"
"       present in the file.\n"
"   --preserve-sybyltypes\n"
"       By default, msys2mol2 will assign sybyl types to each atom using\n"
"       the AssignSybylTypes function.  This option causes msys2mol2 to\n"
"       use the sybyl atom and bond types in the file.\n"
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
        } else if ( !strcmp(pargv[0][j], "--preserve-bondorder")) {
            assign_bondorder = false;
        } else if ( !strcmp(pargv[0][j], "--preserve-sybyltypes")) {
            assign_sybyltypes = false;
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
    SystemPtr h = Load(argv[1]);
    if (assign_bondorder)  AssignBondOrderAndFormalCharge(h);
    if (assign_sybyltypes) AssignSybylTypes(h);
    ExportMol2(h, argv[2], Provenance::fromArgs(argc,argv));
    return 0;
}
