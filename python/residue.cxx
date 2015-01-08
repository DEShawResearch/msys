#include "wrap_obj.hxx"

namespace {
    std::string get_name(residue_t const& r) { return r.name; }
    void set_name(residue_t& r, std::string const& s) { r.name = s; }

    std::string get_insertion(residue_t const& r) { return r.insertion; }
    void set_insertion(residue_t& r, std::string const& s) { r.insertion = s; }
}

namespace desres { namespace msys { 

    void export_residue() {

        class_<residue_t>("residue_t", no_init)
            .add_property("name", get_name, set_name)
            .add_property("insertion", get_insertion, set_insertion)
            .def_readwrite("resid", &residue_t::resid)
            .def_readonly("chain",  &residue_t::chain)
            ;
    }
}}

