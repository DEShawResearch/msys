#include "wrap_obj.hxx"

namespace desres { namespace msys { 

    void export_residue() {

        class_<residue_t>("residue_t", no_init)
            .def_readwrite("name",  &residue_t::name)
            .def_readwrite("num",   &residue_t::num)
            .def_readonly("chain",  &residue_t::chain)
            ;
    }
}}

