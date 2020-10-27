#include "pymod.hxx"
#include <pybind11/stl.h>
#include <msys/system.hxx>

namespace desres { namespace msys { 

    void export_chain(module m) {

        class_<atom_t>(m, "atom_t")
            .def_readwrite("fragid",    &atom_t::fragid)
            .def_readonly("residue",    &atom_t::residue)
            .def_readwrite("x",         &atom_t::x)
            .def_readwrite("y",         &atom_t::y)
            .def_readwrite("z",         &atom_t::z)
            .def_readwrite("charge",    &atom_t::charge)
            .def_readwrite("vx",        &atom_t::vx)
            .def_readwrite("vy",        &atom_t::vy)
            .def_readwrite("vz",        &atom_t::vz)
            .def_readwrite("mass",      &atom_t::mass)
            .def_readwrite("atomic_number", &atom_t::atomic_number)
            .def_readwrite("formal_charge", &atom_t::formal_charge)
            .def_readwrite("aromatic",  &atom_t::aromatic)
            .def_property("name", [](atom_t& a) { return a.name.c_str(); }, [](atom_t& a, std::string const& s) { a.name=s; })
            ;

        class_<bond_t>(m, "bond_t")
            .def_readonly("i",      &bond_t::i)
            .def_readonly("j",      &bond_t::j)
            .def_readwrite("order", &bond_t::order)
            .def("other",           &bond_t::other);
            ;

        class_<residue_t>(m, "residue_t")
            .def_property("name", [](residue_t& a) { return a.name.c_str(); }, [](residue_t& a, std::string const& s) { a.name=s; })
            .def_property("insertion", [](residue_t& a) { return a.insertion.c_str(); }, [](residue_t& a, std::string const& s) { a.insertion=s; })
            .def_readwrite("resid", &residue_t::resid)
            .def_readonly("chain",  &residue_t::chain)
            ;

        class_<chain_t>(m, "chain_t")
            .def_readonly("ct", &chain_t::ct)
            .def_readwrite("name", &chain_t::name)
            .def_readwrite("segid", &chain_t::segid)
            ;

        class_<component_t>(m, "component_t")
            .def_property("name", &component_t::name, &component_t::setName)
            .def("keys", &component_t::keys)
            .def("get",  [](component_t& ct, String const& key) {
                    if (!ct.has(key)) { PyErr_SetString(PyExc_KeyError, key.data()); throw error_already_set(); }
                    return from_value_ref(ct.value(key));
                })
            .def("add",  [](component_t& ct, String const& key, object t) { ct.add(key, as_value_type(t)); })
            .def("set",  [](component_t& ct, String const& key, object val) { to_value_ref(val, ct.value(key)); })
            .def("type", [](component_t& ct, String const& key) -> handle { return ct.has(key) ? from_value_type(ct.type(key)) : none(); })
            .def("remove", &component_t::del)
            ;
    }
}}


