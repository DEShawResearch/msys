#include "dms.hxx"
#include "term_table.hxx"
#include <cstdio>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc<4) {
        fprintf(stderr, "Usage: %s input.dms tablename id1 [id2 [...]]\n",
                argv[0]);
        return 1;
    }
    SystemPtr sys = ImportDMS(argv[1]);
    TermTablePtr table = sys->table(argv[2]);
    IdList ids;
    for (int i=3; i<argc; i++) {
        ids.push_back(atol(argv[i]));
    }
    IdList all = table->findWithAll(ids);
    IdList any = table->findWithAny(ids);
    IdList the = table->findExact(ids);
    printf("found %lu terms with all ids\n", all.size());
    printf("found %lu terms with any ids\n", any.size());
    printf("found %lu terms with the ids\n", the.size());
    return 0;
}

