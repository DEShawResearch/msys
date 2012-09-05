#include "load.hxx"
#include "analyze.hxx"
#include "mol2.hxx"
#include <cstdio>

using namespace desres::msys;

static bool assign_bondorder = true;
static bool assign_sybyltypes = true;

void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "--preserve-bondorder")) {
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
        fprintf(stderr, 
                "Usage: %s [options] input_file output.mol2\n", argv[0]);
        fprintf(stderr, "\nOptions:\n");
        fprintf(stderr, "  --preserve-bondorder    don't recompute bond order\n");
        fprintf(stderr, "  --preserve-sybyltypes   don't reassign sybyl types\n");

        exit(1);
    }
    SystemPtr h = Load(argv[1]);
    if (assign_bondorder)  AssignBondOrderAndFormalCharge(h);
    if (assign_sybyltypes) AssignSybylTypes(h);
    ExportMol2(h, argv[2], Provenance::fromArgs(argc,argv));
    return 0;
}
