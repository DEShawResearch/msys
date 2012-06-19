#include "term_table.hxx"
#include "schema.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {

    SystemPtr sys = System::create();
    for (int i=1; i<argc; i++) {
        TermTablePtr table = AddNonbonded(sys,argv[i], "");
        printf("got table %s with %d atoms, %d term props, %d params\n",
                argv[i], table->atomCount(),
                table->termPropCount(),
                table->params()->propCount());
    }
    return 0;
}
