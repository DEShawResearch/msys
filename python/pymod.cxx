#include <boost/python.hpp>

namespace desres { namespace msys {
    void export_atom();
    void export_bond();
    void export_residue();
    void export_chain();
    void export_param();
    void export_system();
    void export_term();
    void export_vector();
}}

BOOST_PYTHON_MODULE(_msys) {
    desres::msys::export_atom();
    desres::msys::export_bond();
    desres::msys::export_residue();
    desres::msys::export_chain();
    desres::msys::export_param();
    desres::msys::export_system();
    desres::msys::export_term();
    desres::msys::export_vector();
}

