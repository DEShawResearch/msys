#include "wrap_obj.hxx"
#include "annotated_system.hxx"

using namespace desres::msys;

namespace {

    list atom_rings(const AnnotatedSystem& sys, Id atom) {
        MultiIdList rings;
        sys.atomRings(atom, rings);
        list L;
        BOOST_FOREACH(const IdList& ring, rings) {
            list L_in;
            BOOST_FOREACH(Id a, ring)
                L_in.append(a);
            L.append(L_in);
        }
        return L;
    }

    list bond_rings(const AnnotatedSystem& sys, Id bond) {
        MultiIdList rings;
        sys.bondRings(bond, rings);
        list L;
        BOOST_FOREACH(const IdList& ring, rings) {
            list L_in;
            BOOST_FOREACH(Id a, ring)
                L_in.append(a);
            L.append(L_in);
        }
        return L;
    }

    list rings(const AnnotatedSystem& sys) {
        MultiIdList rings;
        sys.rings(rings);
        list L;
        BOOST_FOREACH(const IdList& ring, rings) {
            list L_in;
            BOOST_FOREACH(Id a, ring)
                L_in.append(a);
            L.append(L_in);
        }
        return L;
    }
}

namespace desres { namespace msys {

    void export_annotated_system() {
        class_<AnnotatedSystem, AnnotatedSystemPtr>("AnnotatedSystemPtr", no_init)
            .def("__eq__", list_eq<AnnotatedSystemPtr>)
            .def("__ne__", list_ne<AnnotatedSystemPtr>)
            .def("__hash__", obj_hash<AnnotatedSystemPtr>)
            .def("__init__", make_constructor(&AnnotatedSystem::create))
            .def("system", &AnnotatedSystem::system)
            .def("atomAromatic", &AnnotatedSystem::atomAromatic)
            .def("atomHcount", &AnnotatedSystem::atomHcount)
            .def("atomDegree", &AnnotatedSystem::atomDegree)
            .def("atomValence", &AnnotatedSystem::atomValence)
            .def("atomLonePairs", &AnnotatedSystem::atomLonePairs)
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
            ;
    }
}}
