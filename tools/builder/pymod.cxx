#include "builder.hxx"
#include <boost/python.hpp>

using namespace desres::msys::builder;
using namespace boost::python;

BOOST_PYTHON_MODULE(_builder) {

    class_<defs_t>("Defs", init<>())
        .def("import_charmm_topology", &defs_t::import_charmm_topology)
        ;

    def("build", build);
}

