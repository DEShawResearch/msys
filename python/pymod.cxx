#include <boost/python.hpp>
#include "version.hxx"

namespace desres { namespace msys {
    void export_analyze();
    void export_annotated_system();
    void export_propmap();
    void export_atom();
    void export_bond();
    void export_residue();
    void export_chain();
    void export_param();
    void export_system();
    void export_term();
    void export_override();
    void export_io();
    void export_graph();
    void export_inchi();
    void export_spatial_hash();
}}

BOOST_PYTHON_MODULE(_msys) {
    boost::python::scope().attr("version")=std::string(MSYS_VERSION);
    boost::python::scope().attr("hexversion")=MSYS_VERSION_HEX;
    desres::msys::export_analyze();
    desres::msys::export_annotated_system();
    desres::msys::export_propmap();
    desres::msys::export_atom();
    desres::msys::export_bond();
    desres::msys::export_residue();
    desres::msys::export_chain();
    desres::msys::export_param();
    desres::msys::export_system();
    desres::msys::export_term();
    desres::msys::export_override();
    desres::msys::export_io();
    desres::msys::export_graph();
    desres::msys::export_inchi();
    desres::msys::export_spatial_hash();
}

