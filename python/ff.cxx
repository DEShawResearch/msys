#include "wrap_obj.hxx"
#include "../src/ff.hxx"

namespace {
    using namespace desres::msys;

    list get_es_scale(ff::Rules const& rules) {
        list L;
        for (auto s : rules.es_scale) L.append(object(s));
        return L;
    }
    void set_es_scale(ff::Rules& rules, list L) {
        Id i,n = len(L);
        rules.es_scale.resize(n);
        for (i=0; i<n; i++) rules.es_scale[i] = extract<double>(L[i]);
    }
    list get_lj_scale(ff::Rules const& rules) {
        list L;
        for (auto s : rules.lj_scale) L.append(object(s));
        return L;
    }
    void set_lj_scale(ff::Rules& rules, list L) {
        Id i,n = len(L);
        rules.lj_scale.resize(n);
        for (i=0; i<n; i++) rules.lj_scale[i] = extract<double>(L[i]);
    }

    void build_tuples(ff::Tuples& self, SystemPtr mol, list ids) {
        ff::build(self, mol, ids_from_python(ids));
    }

    void build_component(ff::Forcefield const& self, ff::Component p, SystemPtr mol, ff::Tuples const& t) {
        switch(p) {
            case ff::Component::exclusions:
                ff::build<ff::Component::exclusions>(mol, self, t);
                break;
            default:
                MSYS_FAIL("Unsupported component");
        };
    }

}

BOOST_PYTHON_MODULE(_ff) {
    using namespace desres::msys::ff;

    class_<Rules>("Rules", init<>())
        .def_readwrite("vdw_func", &Rules::vdw_func)
        .def_readwrite("vdw_rule", &Rules::vdw_rule)
        .def_readwrite("exclusions", &Rules::exclusions)
        .add_property("es_scale", get_es_scale, set_es_scale)
        .add_property("lj_scale", get_lj_scale, set_lj_scale)
        ;

    class_<Tuples>("Tuples", init<>())
        .def("build", build_tuples)
        ;
    
    class_<Forcefield>("Forcefield", init<>())
        .def_readwrite("name", &Forcefield::name)
        .def_readonly("rules", &Forcefield::rules)
        .def("build_component", build_component)
        ;

    enum_<Component>("Component")
        .value("exclusions", Component::exclusions)
        ;

}
