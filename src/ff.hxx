#ifndef desres_msys_ff_hxx
#define desres_msys_ff_hxx

#include "system.hxx"

namespace desres { namespace msys { namespace ff {
    // forcefield-specific routines

    enum class Component {
          undefined
        , bonds
        , angles
        , ureybradley
        , propers
        , impropers
        , vdw1
        , vdw2
        , exclusions
        , mass
        , cmap
        , virtuals
    };

    struct Rules {
        //std::vector<std::string> info;
        std::string vdw_func;
        std::string vdw_rule;
        Id exclusions = 0;
        std::vector<double> es_scale;
        std::vector<double> lj_scale;
        //std::vector<Component> components;
        //std::string nbfix;
    };

    struct Forcefield {
        std::string name;
        Rules rules;
        //TemplateMap templates;
        //ParamMap params;
    };

    struct Tuples {
        IdList      atoms;
        MultiIdList bonds;
        MultiIdList angles;
        MultiIdList dihedrals;
    };

    // construct all bonds, angles, dihedrals from connected atoms in fragment
    void build(Tuples& tuples, SystemPtr mol, IdList const& fragment);


    // build terms for component of type p
    template <Component p>
    void build(SystemPtr mol, Forcefield const& ff, Tuples const& tuples);

}}} // namespace

#endif
