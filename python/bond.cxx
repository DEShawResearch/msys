#include "wrap_obj.hxx"

namespace desres { namespace msys { 

    void export_bond() {

        class_<bond_t>("bond_t", no_init)
            .def_readonly("i",      &bond_t::i)
            .def_readonly("j",      &bond_t::j)
            .def_readwrite("order", &bond_t::order)
            .def("other",           &bond_t::other);
            ;

    }
}}

