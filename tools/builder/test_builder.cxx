#include "../../src/dms.hxx"
#include "../../src/clone.hxx"
#include "builder.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<4) {
        fprintf(stderr, "Usage: %s input.dms output.dms 1.top [2.top ...]\n",
                argv[0]);
        return 1;
    }
    SystemPtr dms = ImportDMS(argv[1], true); /* structure only */

    builder::defs_t defs;
    for (int i=3; i<argc; i++) defs.import_charmm_topology(argv[i]);

    for (Id i=0; i<dms->chainCount(); i++) {
        printf("building chain '%s'\n", dms->chain(i).name.c_str());
        builder::build(defs, dms, i);
    }

    dms = Clone(dms, dms->orderedIds());
    ExportDMS(dms, argv[2], Provenance::fromArgs(argc, argv));
    return 0;
}
