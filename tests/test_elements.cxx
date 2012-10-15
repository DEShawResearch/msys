#include "elements.hxx"
#include <stdlib.h>
#include <stdio.h>

using namespace desres::msys;

int main(int argc, char *argv[]) {
    for (int i=0; i<100; i++) {
        printf("%-2s %5.2f\n",
                AbbreviationForElement(i),
                RadiusForElement(i));
    }

    for (int i=1; i<argc; i++) {
        double mass = atof(argv[i]);
        int anum = GuessAtomicNumber(mass);
        printf("mass %f -> anum %d\n", mass, anum);
    }
    return 0;
}

