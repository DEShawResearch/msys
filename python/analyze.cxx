#include <boost/python.hpp>

#include "analyze.hxx"
#include "sssr.hxx"
#include "mol2.hxx"
#include "elements.hxx"
#include "smarts.hxx"

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
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        IdList frags = FindDistinctFragments(mol, fragments);
        list L;
        BOOST_FOREACH(Id frag, frags) L.append(frag);
        return L;
    }

    list find_matches(SmartsPattern const& s, AnnotatedSystemPtr sys, 
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
        /* Yes, we have two interfaces for SSSR, this one and the one in
         * AnnotatedSystem.  This one lets you specify which atoms you
         * want the rings for, and doesn't force you to do any annotation, 
         * which is what we want.  AnnotatedSystem's rings() method only
         * lets you find rings connected to specific atoms or bonds. 
         */
        def("GetSSSR", get_sssr);
        def("ComputeTopologicalIds", ComputeTopologicalIds);
        def("GuessBondConnectivity", GuessBondConnectivity);
        def("FindDistinctFragments", find_distinct_fragments);
        def("RadiusForElement", RadiusForElement);
        def("MassForElement", MassForElement);
        def("PeriodForElement", PeriodForElement);
        def("GroupForElement", GroupForElement);
        def("ElementForAbbreviation", ElementForAbbreviation);
        def("AbbreviationForElement", AbbreviationForElement);
        def("AddHydrogens", AddHydrogens);
        def("GuessHydrogenPositions", GuessHydrogenPositions);
        def("GuessAtomicNumber", GuessAtomicNumber);

        class_<SmartsPattern>("SmartsPattern", init<std::string const&>())
            .def("atomCount", &SmartsPattern::atomCount)
            .def("pattern",   &SmartsPattern::pattern,
                    return_value_policy<copy_const_reference>())
            .def("warnings",  &SmartsPattern::warnings,
                    return_value_policy<copy_const_reference>())
            .def("findMatches",     &find_matches)
            .def("match",     &SmartsPattern::match)
            ;
    }
}}

