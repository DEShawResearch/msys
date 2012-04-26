#include "wrap_obj.hxx"

using namespace desres::msys;

namespace desres { namespace msys { 

    void export_atom() {

        class_<atom_t>("atom_t", no_init)
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
            .def_readwrite("name",      &atom_t::name)
            ;

    }
}}

