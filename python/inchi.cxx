#include "wrap_obj.hxx"
#include "inchi.hxx"

namespace desres { namespace msys { 

    void export_inchi() {

        object cls = class_<InChI>("InChI", no_init)
            .def("create", &InChI::create).staticmethod("create")
            .def("string",  &InChI::string, return_const())
            .def("auxinfo", &InChI::auxinfo, return_const())
            .def("message", &InChI::message, return_const())
            .def("key",     &InChI::key)
            .def("ok",      &InChI::ok)
            ;
        cls.attr("DoNotAddH") = InChI::DoNotAddH;
        cls.attr("SNon") = InChI::SNon;
        cls.attr("FixedH") = InChI::FixedH;
    }
}}

