#include "pymod.hxx"
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <msys/system.hxx>

namespace desres { namespace msys { 

    struct Atom {
        SystemPtr mol;
        Id id;

        bool operator==(Atom const& rhs) const { return mol==rhs.mol && id==rhs.id; }
        bool operator!=(Atom const& rhs) const { return mol!=rhs.mol || id!=rhs.id; }
    };

    void export_chain(module m) {

        class_<Atom>(m, "Atom")
        .def(init<SystemPtr, Id>())
        .def(self == self)
        .def(self != self)
        .def("__hash__", [](Atom const& a) { return hash(make_tuple(a.mol.get(), a.id)); })
        .def("__lt__", [](Atom& lhs, Atom& rhs) { return (lhs.mol < rhs.mol) || (lhs.id < rhs.id); })
        .def_readonly("_ptr", &Atom::mol)
        .def_readonly("id", &Atom::id)
        .def_property("x", [](Atom& a) { return a.mol->atom(a.id).x; },
                [](Atom& a, double val) { a.mol->atom(a.id).x = val; },
                "x component of position")
        .def_property("y", [](Atom& a) { return a.mol->atom(a.id).y; },
                [](Atom& a, double val) { a.mol->atom(a.id).y = val; },
                "y component of position")
        .def_property("z", [](Atom& a) { return a.mol->atom(a.id).z; },
                [](Atom& a, double val) { a.mol->atom(a.id).z = val; },
                "z component of position")
        .def_property("vx", [](Atom& a) { return a.mol->atom(a.id).vx; },
                [](Atom& a, double val) { a.mol->atom(a.id).vx = val; },
                "x component of velocity")
        .def_property("vy", [](Atom& a) { return a.mol->atom(a.id).vy; },
                [](Atom& a, double val) { a.mol->atom(a.id).vy = val; },
                "y component of velocity")
        .def_property("vz", [](Atom& a) { return a.mol->atom(a.id).vz; },
                [](Atom& a, double val) { a.mol->atom(a.id).vz = val; },
                "z component of velocity")
        .def_property("mass", [](Atom& a) { return a.mol->atom(a.id).mass; },
                [](Atom& a, double val) { a.mol->atom(a.id).mass = val; },
                "mass")
        .def_property("charge", [](Atom& a) { return a.mol->atom(a.id).charge; },
                [](Atom& a, double val) { a.mol->atom(a.id).charge= val; },
                "charge")
        .def_property("atomic_number", [](Atom& a) { return a.mol->atom(a.id).atomic_number; },
                [](Atom& a, int val) { a.mol->atom(a.id).atomic_number = val; },
                "atomic number")
        .def_property("formal_charge", [](Atom& a) { return a.mol->atom(a.id).formal_charge; },
                [](Atom& a, int val) { a.mol->atom(a.id).formal_charge = val; },
                "formal charge")
        .def_property("name", [](Atom& a) { return str(a.mol->atom(a.id).name); },
                [](Atom& a, std::string const& val) { a.mol->atom(a.id).name = val; },
                "name")
        .def_property_readonly("fragid", [](Atom& a) { return a.mol->atom(a.id).fragid; },
                "fragment id")
        .def_property_readonly("resid", [](Atom& a) { return a.mol->atom(a.id).residue; },
                "residue id")
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


