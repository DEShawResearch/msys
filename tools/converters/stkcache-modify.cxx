#include <msys/molfile/dtrplugin.hxx>
#include <msys/system.hxx>
#include <fstream>

using namespace desres::msys;
using namespace desres::molfile;

static void usage(FILE *fd) {
    fprintf(fd,
"**stkcache-modify input_v8_stk output_v8_stk_cache input_dtr_name output_dtr_name**\n"
"\n"
"  Modifies a version 8 stk cache file for the input stk, to change the path of\n"
"  the specified dtr file.  Used by desjob and yas.  Not intended to be run\n"
"  manually.\n"
"\n"
);
}

int main(int argc, char **argv) {
    for (int i=0; i< argc; i++) {
        if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
            usage(stdout);
            exit(0);
	}
    }

    if (argc != 5) {
	usage(stderr);
	exit(1);
    }

    char *original_stk = argv[1];
    char *new_stk_cache = argv[2];
    char *original_dtr = argv[3];
    char *new_dtr = argv[4];

    StkReader stk(original_stk);
    stk.init();

    int replaced = 0;
    for (Id i=0, n=stk.nframesets(); i<n; i++) {
        auto dtr = stk.frameset(i);
	if (!strcmp(dtr->path().data(), original_dtr)) {
	    dtr->set_path(new_dtr);
	    replaced++;
	}
    }

    if (replaced == 1) {
	stk.write_cachefile(new_stk_cache);
	exit(0);
    } else if (replaced > 1) {
	fprintf(stderr, "Found multiple instances of %s in stk file %s - cowardly exiting.\n", original_dtr, original_stk);
	exit(1);
    } else {
	fprintf(stderr, "Did not find any instances of %s in stk file %s - cowardly exiting.\n", original_dtr, original_stk);
	exit(1);
    }
}

