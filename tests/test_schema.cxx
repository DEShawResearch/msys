#include "term_table.hxx"
#include "schema.hxx"
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {

    SystemPtr sys = System::create();
    sys->addChain();
    sys->addResidue(0);
    for (int i=1; i<argc; i++) {
        TermTablePtr table = AddTable(sys,argv[i], "");
        printf("got table %s with %d atoms, %d term props, %d params\n",
                argv[i], table->atomCount(),
                table->termPropCount(),
                table->params()->propCount());
        IdList ids(table->atomCount());
        for (Id j=0; j<table->atomCount(); j++) {
            ids[j] = sys->addAtom(0);
        }
        Id term = table->addTerm(ids,BadId);
        /* no such param property */
        try {
            table->propValue(term, "x0");
        } catch (std::exception& e) {
            printf("message for bad name in propValue: %s\n", e.what());
        }
        /* no such term for param */
        try {
            table->propValue(32, "fcx");
        } catch (std::exception& e) {
            printf("message for bad term in propValue: %s\n", e.what());
        }
        /* no such term property */
        try {
            table->termPropValue(term, "fcx");
        } catch (std::exception& e) {
            printf("message for bad name in termPropValue: %s\n", e.what());
        }
        /* no such term for term prop */
        try {
            table->termPropValue(32, "x0");
        } catch (std::exception& e) {
            printf("message for bad term in termPropValue: %s\n", e.what());
        }
    }
    return 0;
}
