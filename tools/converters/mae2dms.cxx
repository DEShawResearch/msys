#include "mae.hxx"
#include <dms.hxx>
#include <boost/filesystem.hpp>

#include <string.h>
#include <fstream>

using namespace desres::msys;

static void write_provenance( int argc, char * argv[], SystemPtr h ) {
    std::string path = boost::filesystem::system_complete(".").string();
    std::string cmdline;
    for (int i=0; i<argc; i++) {
        cmdline += argv[i];
        if (i!=argc-1) cmdline += " ";
    }
    ParamTablePtr p = ParamTable::create();
    p->addProp("workdir", StringType);
    p->addProp("cmdline", StringType);
    Id param = p->addParam();
    p->value(param,0)=path;
    p->value(param,1)=cmdline;
    h->addExtra("provenance", p);
}

static bool ignore_unrecognized = false;

void parse_cmdline( int *pargc, char ***pargv ) {
    int i,j=0;
                           
    for (i=0; i<pargc[0]; i++) {
        pargv[0][j] = pargv[0][i];
        if ( !strcmp(pargv[0][j], "--ignore-unrecognized")) {
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
        fprintf(stderr, "\nUsage: %s input.mae output.dms [options]\n", 
                argv[0]);
        fprintf(stderr, "\nOptions:\n\n");
        fprintf(stderr, "  --ignore-unrecognized   -- skip unrecognized ffio_ff subblocks\n");
        fprintf(stderr, "\n");
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

    write_provenance(argc, argv, h);

    try {
        ExportDMS(h, argv[2]);
    }
    catch (std::exception& e) {
        fprintf(stderr, "mae2dms: Failed exporting dms to '%s'\n", argv[2]);
        fprintf(stderr, "ERROR: %s\n", e.what());
        exit(1);
    }
    return 0;
}

