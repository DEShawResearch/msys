#include "wrap_obj.hxx"
#include "annotated_system.hxx"
#include "sssr.hxx"

using namespace desres::msys;

namespace {

    list atom_rings(const AnnotatedSystem& sys, Id atom) {
        MultiIdList rings;
        sys.atomRings(atom, rings);
        return to_python(rings);
    }

    list bond_rings(const AnnotatedSystem& sys, Id bond) {
        MultiIdList rings;
        sys.bondRings(bond, rings);
        return to_python(rings);
    }

    list rings(const AnnotatedSystem& sys) {
        MultiIdList rings;
        sys.rings(rings);
        return to_python(rings);
    }
    list errors(AnnotatedSystem const& a) {
        list L;
        for (auto s : a.errors()) {
            L.append(object(s));
        }
        return L;
    }

    list wrap_atoms(AnnotatedSystem const& a) {
        return to_python(a.atoms());
    }

}

namespace desres { namespace msys {

    void export_annotated_system() {

        enum_<AnnotatedSystem::Flags>("AnnotatedSystemFlags")
            .value("Default",           AnnotatedSystem::Default)
            .value("AllowBadCharges",   AnnotatedSystem::AllowBadCharges)
            ;

        class_<AnnotatedSystem, boost::noncopyable>("AnnotatedSystem", init<SystemPtr, unsigned>())
            .def("atoms", wrap_atoms)
            .def("atomAromatic", &AnnotatedSystem::atomAromatic)
            .def("atomHcount", &AnnotatedSystem::atomHcount)
            .def("atomDegree", &AnnotatedSystem::atomDegree)
            .def("atomValence", &AnnotatedSystem::atomValence)
            .def("atomLoneElectrons", &AnnotatedSystem::atomLoneElectrons)
            .def("atomHybridization", &AnnotatedSystem::atomHybridization)
            .def("atomRingBonds", &AnnotatedSystem::atomRingBonds)
            .def("atomRingCount", &AnnotatedSystem::atomRingCount)
            .def("atomInRingSize", &AnnotatedSystem::atomInRingSize)
            .def("atomRings", atom_rings)
            .def("bondAromatic", &AnnotatedSystem::bondAromatic)
            .def("bondRingCount", &AnnotatedSystem::bondRingCount)
            .def("bondInRingSize", &AnnotatedSystem::bondInRingSize)
            .def("bondRings", bond_rings)
            .def("ringCount", &AnnotatedSystem::ringCount)
            .def("rings", rings)
            .def("errors", errors)
            ;
    }
}}
