#include "wrap_obj.hxx"
#include "inchi.hxx"

namespace desres { namespace msys { 

    void export_inchi() {

#ifndef MSYS_WITHOUT_INCHI
        scope cls = class_<InChI>("InChI", no_init)
            .def("create", &InChI::create).staticmethod("create")
            .def("string",  &InChI::string, return_const())
            .def("auxinfo", &InChI::auxinfo, return_const())
            .def("message", &InChI::message, return_const())
            .def("key",     &InChI::key)
            .def("ok",      &InChI::ok)
            ;

        // since we assigned cls to a scope above, this enum is in class scope
        enum_<InChI::Flags>("Flags")
            .value("DoNotAddH", InChI::DoNotAddH)
            .value("SNon",      InChI::SNon)
            .value("FixedH",    InChI::FixedH)
            ;
#endif
    }
}}

