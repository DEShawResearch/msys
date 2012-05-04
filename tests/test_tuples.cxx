#include "tuples.hxx"
#include "dms.hxx"

using namespace desres::msys;

int main(int argc, char *argv[]) {
    if (argc>1) {
        SystemPtr sys = ImportDMS(argv[1]);
        TermTablePtr nonbonded = sys->table("nonbonded");
        ParamTablePtr combined = sys->auxTable("nonbonded_combined_param");
        TermTablePtr tuples = sys->addTable("nonbonded_combined", 2, 
                nonbonded->params());

        CreateTuplesFromCombined(nonbonded, combined, tuples);
    }
    return 0;
}

