#include "schema.hxx"
#include "../fastjson/parse.hxx"
#include "../fastjson/print.hxx"

using namespace desres::msys;
using desres::msys::fastjson::Json;

static char schemas_as_json[] = R"SCHEMA(
[
     [ "angle_harm",   "bond", 3, ["theta0", "fc"], [["constrained", 0]]],
     [ "angle_fbhw",   "bond", 3, ["sigma", "theta0", "fc"]],
     [ "dihedral_trig", "bond",4, ["phi0","fc0","fc1","fc2","fc3","fc4","fc5","fc6"]],
     [ "dihedral6_trig","bond",6, ["phi0","fc0","fc2","fc4"]],
     [ "improper_anharm", "bond",4, ["fc2","fc4"]],
     [ "improper_fbhw", "bond",4, ["sigma","phi0","fc"]],
     [ "improper_harm", "bond",4, ["phi0","fc"]],
     [ "inplanewag_harm","bond",4,["w0","fc"]],
     [ "pair_12_6_es", "bond", 2, ["aij", "bij", "qij"]],
     [ "pair_exp_6_es","bond", 2, ["aij", "bij", "cij", "qij"]],
     [ "posre_harm",   "bond", 1, ["fcx", "fcy", "fcz"], ["x0", "y0", "z0"]],
     [ "posre_fbhw",   "bond", 1, ["fc", "sigma"], ["x0", "y0", "z0"]],
     [ "stretch_harm", "bond", 2, ["r0", "fc"],     [["constrained", 0]]],
     [ "stretch_morse", "bond", 2,["r0", "d", "a"]],
     [ "softstretch_harm", "bond", 2, ["r0", "fc", "alpha"]],
     [ "torsiontorsion_cmap", "bond", 8, [["cmapid", 2]]],
     [ "pseudopol_fermi","bond",4,["a", "b", "cutoff"]],

     [ "alchemical_stretch_harm", "bond", 2, ["r0A", "fcA", "r0B", "fcB"], [["schedule", 2]]],
     [ "alchemical_stretch_morse", "bond", 2,["r0A", "dA", "aA", "r0B", "dB", "aB"], [["schedule", 2]]],
     [ "alchemical_softstretch_harm", "bond", 2, ["r0A", "fcA", "alphaA", "r0B", "fcB", "alphaB"], [["schedule", 2]]],
     [ "alchemical_angle_harm", "bond", 3, ["theta0A", "fcA", "theta0B", "fcB"], [["schedule", 2]]],
     [ "alchemical_angle_harm_soft", "bond", 3, ["theta0A", "fcA", "epsilonA", "theta0B", "fcB", "epsilonB"], [["schedule", 2]]],
     [ "alchemical_dihedral_trig", "bond",4, 
        ["phi0A","fc0A","fc1A","fc2A","fc3A","fc4A","fc5A","fc6A",
         "phi0B","fc0B","fc1B","fc2B","fc3B","fc4B","fc5B","fc6B"], [["schedule", 2]]],
     [ "alchemical_dihedral_trig_soft", "bond",4, 
        ["phi0A","fc0A","fc1A","fc2A","fc3A","fc4A","fc5A","fc6A","epsilonA",
         "phi0B","fc0B","fc1B","fc2B","fc3B","fc4B","fc5B","fc6B","epsilonB"], [["schedule", 2]]],
     [ "alchemical_pair_12_6_es", "bond", 2, 
        ["aijA", "bijA", "qijA",
         "aijB", "bijB", "qijB"], [["schedule", 2]]],
     [ "alchemical_pair_exp_6_es","bond", 2, 
        ["aijA", "bijA", "cijA", "qijA",
         "aijB", "bijB", "cijB", "qijB"], [["schedule", 2]]],
     [ "alchemical_improper_harm","bond", 4,
        ["phi0A", "fcA", "phi0B", "fcB"], [["schedule", 2]]],
     [ "alchemical_improper_harm_soft","bond", 4,
        ["phi0A", "fcA", "epsilonA", "phi0B", "fcB", "epsilonB"], [["schedule", 2]]],
     [ "alchemical_torsiontorsion_cmap", "bond", 8, [["cmapidA", 2], ["cmapidB", 2]], [["schedule", 2]]],

     [ "constraint_hoh", "constraint", 3, ["theta","r1","r2"]],
     [ "constraint_ah1", "constraint", 2, ["r1"]],
     [ "constraint_ah2", "constraint", 3, ["r1","r2"]],
     [ "constraint_ah3", "constraint", 4, ["r1","r2","r3"]],
     [ "constraint_ah4", "constraint", 5, ["r1","r2","r3","r4"]],
     [ "constraint_ah5", "constraint", 6, ["r1","r2","r3","r4","r5"]],
     [ "constraint_ah6", "constraint", 7, ["r1","r2","r3","r4","r5","r6"]],
     [ "constraint_ah7", "constraint", 8, ["r1","r2","r3","r4","r5","r6","r7"]],
     [ "constraint_ah8", "constraint", 9, ["r1","r2","r3","r4","r5","r6","r7","r8"]],

     [ "rigid_explicit2", "constraint", 2, ["x0", "y0", "z0", "x1", "y1", "z1"]],
     [ "rigid_explicit3", "constraint", 3, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2"]],
     [ "rigid_explicit4", "constraint", 4, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3"]],
     [ "rigid_explicit5", "constraint", 5, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4"]],
     [ "rigid_explicit6", "constraint", 6, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "x5", "y5", "z5"]],
     [ "rigid_explicit7", "constraint", 7, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "x5", "y5", "z5", "x6", "y6", "z6"]],
     [ "rigid_explicit8", "constraint", 8, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "x5", "y5", "z5", "x6", "y6", "z6", "x7", "y7", "z7"]],
     [ "rigid_explicit9", "constraint", 9, ["x0", "y0", "z0", "x1", "y1", "z1", "x2", "y2", "z2", "x3", "y3", "z3", "x4", "y4", "z4", "x5", "y5", "z5", "x6", "y6", "z6", "x7", "y7", "z7", "x8", "y8", "z8"]],

     [ "virtual_fdat3",   "virtual", 4, [["c1",1,"distance"],
                                      ["c2",1,"angle"],
                                      ["c3",1,"torsion"]]],

     [ "virtual_lc2",     "virtual", 3, ["c1"]],
     [ "virtual_lc2n",    "virtual", 3, ["c1"]],
     [ "virtual_lc3",     "virtual", 4, ["c1","c2"]],
     [ "virtual_lc3n",    "virtual", 4, ["c1","c2"]],
     [ "virtual_lc4",     "virtual", 5, ["c1","c2","c3"]],
     [ "virtual_lc4n",    "virtual", 5, ["c1","c2","c3"]],
     [ "virtual_lc5",     "virtual", 6, ["c1","c2","c3","c4"]],
     [ "virtual_lc5n",    "virtual", 6, ["c1","c2","c3","c4"]],
     [ "virtual_lc6",     "virtual", 7, ["c1","c2","c3","c4","c5"]],
     [ "virtual_lc6n",    "virtual", 7, ["c1","c2","c3","c4","c5"]],
     [ "virtual_lc7",     "virtual", 8, ["c1","c2","c3","c4","c5","c6"]],
     [ "virtual_lc7n",    "virtual", 8, ["c1","c2","c3","c4","c5","c6"]],
     [ "virtual_midpoint","virtual", 3, ["c1"]], /* identical to lc2 */
     [ "virtual_out3",    "virtual", 4, ["c1","c2","c3"]],
     [ "virtual_out3n",   "virtual", 4, ["c1","c2","c3"]],
     [ "virtual_sp3",     "virtual", 4, ["c1","c2"]],

     [ "exclusion", "exclusion", 2 ]
])SCHEMA";

static char schemas_nb_as_json[] = R"SCHEMA(
[
     [ "vdw_12_6", "nonbonded", 1, ["sigma", "epsilon"]],
     [ "vdw_exp_6", "nonbonded", 1, ["alpha", "epsilon", "rmin"]], 
     [ "vdw_exp_6s", "nonbonded",1, ["sigma", "epsilon", "lne"]],
     [ "polynomial_cij", "nonbonded", 1, 
        [ "c1",  "c2",  "c3",  "c4",  "c5",  "c6",  "c7",  "c8", 
          "c9",  "c10", "c11", "c12", "c13", "c14", "c15", "c16"]] 
])SCHEMA";

namespace {
    Json schemas, schemas_nb;
    struct _ {
        _() {
            std::istringstream in(schemas_as_json);
            desres::msys::fastjson::parse_json(in, schemas);
            std::istringstream nb(schemas_nb_as_json);
            desres::msys::fastjson::parse_json(nb, schemas_nb);
        }
    } initializer;

    void parse_column(Json const& js, column_t* col) {
        if (js.kind()==Json::String) {
            col->name = js.as_string();
        } else {
            col->name = js.elem(0).as_string();
            if (js.size()>1) {
                col->type = js.elem(1).as_int();
                if (js.size()>2) {
                    col->comment = js.elem(2).as_string();
                }
            }
        }
    }
    schema_t make_schema(Json const& js) {
        schema_t s;
        s.name = js.elem(0).as_string();
        s.category = js.elem(1).as_string();
        s.nsites = js.elem(2).as_int();
        if (js.size()>3) {
            Json const& params = js.elem(3);
            for (int i=0; i<params.size(); i++) {
                parse_column(params.elem(i), s.param_props+i);
            }
        }
        if (js.size()>4) {
            Json const& term_props = js.elem(4);
            for (int i=0; i<term_props.size(); i++) {
                parse_column(term_props.elem(i), s.term_props+i);
            }
        }
        return s;
    }
}

namespace desres { namespace msys { 

    unsigned schema_count() { return schemas.size(); }
    unsigned nonbonded_schema_count() { return schemas_nb.size(); }

    schema_t schema(unsigned i) {
        return make_schema(schemas.elem(i));
    }
    schema_t nonbonded_schema(unsigned i) {
        return make_schema(schemas_nb.elem(i));
    }

}}

