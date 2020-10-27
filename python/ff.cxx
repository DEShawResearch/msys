#include "pymod.hxx"
#include <pybind11/stl.h>
#include <msys/ff.hxx>

using namespace desres::msys::ff;

namespace desres { namespace msys {

void export_ff(module m) {

    class_<Rules>(m, "Rules")
        .def(init<>())
        .def_readwrite("vdw_func", &Rules::vdw_func)
        .def_readwrite("vdw_rule", &Rules::vdw_rule)
        .def_readwrite("exclusions", &Rules::exclusions)
        .def_readwrite("es_scale", &Rules::es_scale)
        .def_readwrite("lj_scale", &Rules::lj_scale)
        ;

    class_<Tuples>(m, "Tuples")
        .def(init<>())
        .def("build", [](Tuples& self, SystemPtr mol, IdList const& fragment) { build(self, mol, fragment); })
        ;
    
    class_<Forcefield>(m, "Forcefield")
        .def(init<>())
        .def_readwrite("name", &Forcefield::name)
        .def_readonly("rules", &Forcefield::rules)
        .def("build_component", [](Forcefield const& self, Component p, SystemPtr mol, Tuples const& t) {
            switch(p) {
                case Component::exclusions:
                    build<Component::exclusions>(mol, self, t);
                    break;
                default:
                    MSYS_FAIL("Unsupported component");
            };
            })
        ;

    enum_<Component>(m, "Component")
        .value("exclusions", Component::exclusions)
        ;

}

}}
