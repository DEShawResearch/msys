#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "inchi.hxx"

using namespace pybind11;

namespace desres { namespace msys { 

    void export_inchi(module m) {

#ifndef MSYS_WITHOUT_INCHI
        auto cls = class_<InChI>(m, "InChI")
            .def_static("create", &InChI::create)
            .def("string",  &InChI::string)
            .def("auxinfo", &InChI::auxinfo)
            .def("message", &InChI::message)
            .def("key",     &InChI::key)
            .def("ok",      &InChI::ok)
            ;

        enum_<InChI::Flags>(cls, "Flags", arithmetic())
            .value("DoNotAddH", InChI::DoNotAddH)
            .value("SNon",      InChI::SNon)
            .value("FixedH",    InChI::FixedH)
            ;
#endif
    }
}}

