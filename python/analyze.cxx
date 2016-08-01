#include <boost/python.hpp>

#include "wrap_obj.hxx"
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

    list find_distinct_fragments(SystemPtr mol) {
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        return to_python(FindDistinctFragments(mol, fragments));
    }

    list find_matches(SmartsPattern const& s, AnnotatedSystem const& sys, 
                                              list const& starts) {
        return to_python(s.findMatches(sys, ids_from_python(starts)));
    }

    list wrap_sssr(SystemPtr mol, IdList const& atoms, bool all_relevant=false)
    {
        return to_python(GetSSSR(mol, atoms, all_relevant));
    }
    list ring_systems(SystemPtr mol, IdList const& atoms) {
        return to_python(RingSystems(mol, GetSSSR(mol, atoms, true)));
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
        def("GetSSSR", wrap_sssr);
        def("RingSystems", ring_systems);
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
        def("Analyze", Analyze);
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

