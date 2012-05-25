#include "schema.hxx"

using namespace desres::msys;

static schema_t schemas[] = {
    /* bond terms */
    { "angle_harm",   "bond", 3, {"theta0", "fc"}, {{"constrained", 0}}},
    { "dihedral_trig", "bond",4, {"phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"}},
    { "dihedral6_trig","bond",6, {"phi0","fc0","fc2","fc4"}},
    { "improper_anharm", "bond",4, {"fc2","fc4"}},
    { "improper_harm", "bond",4, {"phi0","fc"}},
    { "inplanewag_harm","bond",4,{"w0","fc"}},
    { "pair_12_6_es", "bond", 2, {"aij", "bij", "qij"}},
    { "pair_exp_6_es","bond", 2, {"aij", "bij", "cij", "qij"}},
    { "posre_harm",   "bond", 1, {"fcx", "fcy", "fcz"}, {"x0", "y0", "z0"}},
    { "stretch_harm", "bond", 2, {"r0", "fc"},     {{"constrained", 0}}},
    { "torsiontorsion_cmap", "bond", 8, {{"cmapid", 2}}},
    { "pseudopol_fermi","bond",4,{"a", "b", "cutoff"}},

    /* alchemical bond terms */
    { "alchemical_stretch_harm", "bond", 2, {"r0A", "fcA", "r0B", "fcB"}},
    { "alchemical_angle_harm", "bond", 3, {"theta0A", "fcA", "theta0B", "fcB"}},
    { "alchemical_dihedral_trig", "bond",4, 
        {"phi0A","fc0A","fc1A","fc2A","fc3A","fc4A","fc5A","fc6A",
         "phi0B","fc0B","fc1B","fc2B","fc3B","fc4B","fc5B","fc6B"}},
    { "alchemical_pair_12_6_es", "bond", 2, 
        {"aijA", "bijA", "qijA",
         "aijB", "bijB", "qijB"}},
    { "alchemical_pair_exp_6_es","bond", 2, 
        {"aijA", "bijA", "cijA", "qijA",
         "aijB", "bijB", "cijB", "qijB"}},
    { "alchemical_improper_harm","bond", 4,
        {"phi0A", "fcA", "phi0B", "fcB"}},

    /* constraints */
    { "constraint_hoh", "constraint", 3, {"theta","r1","r2"}},
    { "constraint_ah1", "constraint", 2, {"r1"}},
    { "constraint_ah2", "constraint", 3, {"r1","r2"}},
    { "constraint_ah3", "constraint", 4, {"r1","r2","r3"}},
    { "constraint_ah4", "constraint", 5, {"r1","r2","r3","r4"}},
    { "constraint_ah5", "constraint", 6, {"r1","r2","r3","r4","r5"}},
    { "constraint_ah6", "constraint", 7, {"r1","r2","r3","r4","r5","r6"}},
    { "constraint_ah7", "constraint", 8, {"r1","r2","r3","r4","r5","r6","r7"}},
    { "constraint_ah8", "constraint", 9, {"r1","r2","r3","r4","r5","r6","r7","r8"}},

    /* virtuals */
    { "virtual_fdat3",   "virtual", 4, {{"c1",1,"distance"},
                                      {"c2",1,"angle"},
                                      {"c3",1,"torsion"}}},

    { "virtual_lc2",     "virtual", 3, {"c1"}},
    { "virtual_lc2n",    "virtual", 3, {"c1"}},
    { "virtual_lc3",     "virtual", 4, {"c1","c2"}},
    { "virtual_lc3n",    "virtual", 4, {"c1","c2"}},
    { "virtual_lc4",     "virtual", 5, {"c1","c2","c3"}},
    { "virtual_lc4n",    "virtual", 5, {"c1","c2","c3"}},
    { "virtual_midpoint","virtual", 3, {"c1"}}, /* identical to lc2 */
    { "virtual_out3",    "virtual", 4, {"c1","c2","c3"}},
    { "virtual_out3n",   "virtual", 4, {"c1","c2","c3"}},
    { "virtual_sp3",     "virtual", 4, {"c1","c2"}},

    /* exclusion */
    { "exclusion", "exclusion", 2 }
};

static schema_t nonbonded_schemas[] = {
    { "vdw_12_6", "nonbonded", 1, {"sigma", "epsilon"}},
        //{"geometric", "arithmetic", "arithmetic/geometric" }},
    { "vdw_exp_6", "nonbonded", 1, {"alpha", "epsilon", "rmin"}}, 
        //{"lb/geometric"}},
    { "vdw_exp_6s", "nonbonded",1, {"sigma", "epsilon", "lne"}},
        //{"lb/geometric"}},
    { "polynomial_cij", "nonbonded", 1, 
        { "c1",  "c2",  "c3",  "c4",  "c5",  "c6",  "c7",  "c8", 
          "c9",  "c10", "c11", "c12", "c13", "c14", "c15", "c16"}} 
        //{"geometric", "arithmetic", "arithmetic/geometric" }},
};

namespace desres { namespace msys { 

    unsigned schema_count() { 
        return sizeof(schemas)/sizeof(schemas[0]);
    }
    const char* schema_name(unsigned i) {
        return schemas[i].name;
    }
    unsigned nonbonded_schema_count() {
        return sizeof(nonbonded_schemas)/sizeof(nonbonded_schemas[0]);
    }
    const char* nonbonded_schema_name(unsigned i) {
        return nonbonded_schemas[i].name;
    }

    unsigned nonbonded_schema_count();

    const schema_t* find_schema(const std::string& name) {
        unsigned i,n = sizeof(schemas)/sizeof(schemas[0]);
        for (i=0; i<n; i++) {
            if (schemas[i].name && name==schemas[i].name) return schemas+i;
        }
        return NULL;
    }
    const schema_t* find_nonbonded_schema(const std::string& name) {
        unsigned i,n = sizeof(nonbonded_schemas)/sizeof(nonbonded_schemas[0]);
        for (i=0; i<n; i++) {
            if (nonbonded_schemas[i].name && name==nonbonded_schemas[i].name) 
                return nonbonded_schemas+i;
        }
        return NULL;
    }


}}

