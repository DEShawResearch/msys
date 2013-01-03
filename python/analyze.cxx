#include "analyze.hxx"
#include "sssr.hxx"
#include "mol2.hxx"
#include "elements.hxx"
#include "smarts.hxx"
#include <boost/python.hpp>

using namespace desres::msys;
using namespace boost::python;

namespace {
    void assign_1(SystemPtr mol) {
        AssignBondOrderAndFormalCharge(mol);
    }
    void assign_2(SystemPtr mol, IdList const& ids) {
        AssignBondOrderAndFormalCharge(mol, ids);
    }
    void assign_3(SystemPtr mol, IdList const& ids, int total_charge) {
        AssignBondOrderAndFormalCharge(mol, ids, total_charge);
    }
    list get_sssr(SystemPtr mol, IdList const& ids, bool all_relevant) {
        MultiIdList rings = GetSSSR(mol, ids, all_relevant);
        list L;
        for (unsigned i=0; i<rings.size(); i++) L.append(rings[i]);
        return L;
    }
    list find_distinct_fragments(SystemPtr mol) {
        IdList frags = FindDistinctFragments(mol);
        list L;
        BOOST_FOREACH(Id frag, frags) L.append(frag);
        return L;
    }

    list find_matches(SmartsPattern const& s, SystemPtr sys, 
                                              IdList const& starts)
    {
        MultiIdList results = s.findMatches(sys, starts);
        list L;
        BOOST_FOREACH(IdList const& ids, results) {
            list m;
            BOOST_FOREACH(Id id, ids) m.append(id);
            L.append(m);
        }
        return L;
    }

}

namespace desres { namespace msys { 

    void export_analyze() {
        def("AssignBondOrderAndFormalCharge", assign_1);
        def("AssignBondOrderAndFormalCharge", assign_2);
        def("AssignBondOrderAndFormalCharge", assign_3);
        def("AssignSybylTypes", AssignSybylTypes);
        def("GetSSSR", get_sssr);
        def("GuessBondConnectivity", GuessBondConnectivity);
        def("FindDistinctFragments", find_distinct_fragments);
        def("RadiusForElement", RadiusForElement);
        def("MassForElement", MassForElement);

        class_<SmartsPattern>("SmartsPattern", init<std::string const&>())
            .def("atomCount", &SmartsPattern::atomCount)
            .def("pattern",   &SmartsPattern::pattern,
                    return_value_policy<copy_const_reference>())
            .def("warnings",  &SmartsPattern::warnings,
                    return_value_policy<copy_const_reference>())
            .def("findMatches",     &find_matches)
            ;
    }
}}

