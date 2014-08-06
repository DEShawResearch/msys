#include "dms.hxx"
#include "term_table.hxx"
#include <cstdio>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<5) {
        fprintf(stderr, "Usage: %s input.dms tablename property value\n",
                argv[0]);
        return 1;
    }
    SystemPtr sys = ImportDMS(argv[1]);
    TermTablePtr table = sys->table(argv[2]);
    ParamTablePtr params = table->params();
    Id col = params->propIndex(argv[3]);
    for (int i=4; i<argc; i++) {
        IdList ids;
        switch (params->propType(col)) {
            case IntType: ids=params->findInt(col, atoi(argv[i])); break;
            case FloatType: ids=params->findFloat(col, atof(argv[i])); break;
            default:
            case StringType: ids=params->findString(col, argv[i]); break;
        }
        printf("found %4lu params with %s value %s\n", 
                ids.size(), argv[3], argv[i]);
    }
    return 0;
}

