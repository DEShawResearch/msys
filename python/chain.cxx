#include "wrap_obj.hxx"

namespace desres { namespace msys { 

    void export_chain() {

        class_<chain_t>("chain_t", no_init)
            .def_readwrite("name", &chain_t::name)
            ;
    }
}}

