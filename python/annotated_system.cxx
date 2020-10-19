#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "annotated_system.hxx"

using namespace desres::msys;
using namespace pybind11;

namespace desres { namespace msys {

    void export_annotated_system(module m) {

        enum_<AnnotatedSystem::Flags>(m, "AnnotatedSystemFlags")
            .value("Default",           AnnotatedSystem::Default)
            .value("AllowBadCharges",   AnnotatedSystem::AllowBadCharges)
            ;

        class_<AnnotatedSystem>(m, "AnnotatedSystem")
            .def(init<SystemPtr, unsigned>())
            .def("atoms", &AnnotatedSystem::atoms)
            .def("atomAromatic", &AnnotatedSystem::atomAromatic)
            .def("atomHcount", &AnnotatedSystem::atomHcount)
            .def("atomDegree", &AnnotatedSystem::atomDegree)
            .def("atomValence", &AnnotatedSystem::atomValence)
            .def("atomLoneElectrons", &AnnotatedSystem::atomLoneElectrons)
            .def("atomHybridization", &AnnotatedSystem::atomHybridization)
            .def("atomRingBonds", &AnnotatedSystem::atomRingBonds)
            .def("atomRingCount", &AnnotatedSystem::atomRingCount)
            .def("atomInRingSize", &AnnotatedSystem::atomInRingSize)
            .def("atomRings", &AnnotatedSystem::atomRings)
            .def("bondAromatic", &AnnotatedSystem::bondAromatic)
            .def("bondRingCount", &AnnotatedSystem::bondRingCount)
            .def("bondInRingSize", &AnnotatedSystem::bondInRingSize)
            .def("bondRings", &AnnotatedSystem::bondRings)
            .def("ringCount", &AnnotatedSystem::ringCount)
            .def("rings", &AnnotatedSystem::rings)
            .def("errors", &AnnotatedSystem::errors)
            ;
    }
}}
