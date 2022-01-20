#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "analyze.hxx"
#include "sssr.hxx"
#include "mol2.hxx"
#include "elements.hxx"
#include "smarts.hxx"

using namespace desres::msys;
using namespace pybind11;

namespace {
    void assign_1(SystemPtr mol, bool compute_resonant_charges, int timeout_ms) {
        unsigned flags = 0;
        if (compute_resonant_charges) flags |= AssignBondOrder::ComputeResonantCharges;
        AssignBondOrderAndFormalCharge(mol, flags, std::chrono::milliseconds(timeout_ms));
    }
    void assign_2(SystemPtr mol, IdList const& ids, bool compute_resonant_charges, int timeout_ms) {
        unsigned flags = 0;
        if (compute_resonant_charges) flags |= AssignBondOrder::ComputeResonantCharges;
        auto topids = ComputeTopologicalIds(mol);
        AssignBondOrderAndFormalCharge(mol, ids, INT_MAX, flags, topids, std::chrono::milliseconds(timeout_ms));
    }
    void assign_3(SystemPtr mol, IdList const& ids, int total_charge, bool compute_resonant_charges, int timeout_ms) {
        unsigned flags = 0;
        if (compute_resonant_charges) flags |= AssignBondOrder::ComputeResonantCharges;
        auto topids = ComputeTopologicalIds(mol);
        AssignBondOrderAndFormalCharge(mol, ids, total_charge, flags, topids, std::chrono::milliseconds(timeout_ms));
    }


    KekuleResult kekule_1(SystemPtr mol, int timeout_ms) {
        return KekuleStructures(mol, INT_MAX, std::chrono::milliseconds(timeout_ms));
    }
    KekuleResult kekule_2(SystemPtr mol, int total_charge, int timeout_ms) {
        return KekuleStructures(mol, total_charge, std::chrono::milliseconds(timeout_ms));
    }

    std::map<Id,IdList> find_distinct_fragments(SystemPtr mol, std::vector<std::string> const& keys) {
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        return FindDistinctFragments(mol, fragments, keys);
    }

    MultiIdList ring_systems(SystemPtr mol, IdList const& atoms) {
        return RingSystems(mol, GetSSSR(mol, atoms, true));
    }
    double elec_for_element(int n) {
        return DataForElement(n).eneg;
    }

    dict get_bonds_angles_dihedrals(SystemPtr mol) {
        std::vector<IdList> non_pseudo_bonds;
        std::vector<IdList> pseudo_bonds;
        std::vector<IdList> angles;
        std::vector<IdList> dihedrals;
        GetBondsAnglesDihedrals(mol, mol->atoms(), non_pseudo_bonds, pseudo_bonds, angles, dihedrals);
        dict d;
        d["non_pseudo_bonds"] = cast(non_pseudo_bonds);
        d["pseudo_bonds"] = cast(pseudo_bonds);
        d["angles"] = cast(angles);
        d["dihedrals"] = cast(dihedrals);
        return d;
    }

}

namespace desres { namespace msys { 

    void export_analyze(module m) {
        m.def("AssignBondOrderAndFormalCharge", assign_1);
        m.def("AssignBondOrderAndFormalCharge", assign_2);
        m.def("AssignBondOrderAndFormalCharge", assign_3);
        m.def("KekuleStructures", kekule_1);
        m.def("KekuleStructures", kekule_2);
        /* Yes, we have two interfaces for SSSR, this one and the one in
         * AnnotatedSystem.  This one lets you specify which atoms you
         * want the rings for, and doesn't force you to do any annotation, 
         * which is what we want.  AnnotatedSystem's rings() method only
         * lets you find rings connected to specific atoms or bonds. 
         */
        m.def("GetSSSR", GetSSSR);
        m.def("RingSystems", ring_systems);
        m.def("ComputeTopologicalIds", ComputeTopologicalIds);
        m.def("GuessBondConnectivity", GuessBondConnectivity);
        m.def("FindDistinctFragments", find_distinct_fragments);
        m.def("RadiusForElement", RadiusForElement);
        m.def("MassForElement", MassForElement);
        m.def("PeriodForElement", PeriodForElement);
        m.def("GroupForElement", GroupForElement);
        m.def("ElementForAbbreviation", ElementForAbbreviation);
        m.def("GuessHydrogenPositions", GuessHydrogenPositions);
        m.def("AbbreviationForElement", AbbreviationForElement);
        m.def("Analyze", Analyze);
        m.def("GuessAtomicNumber", GuessAtomicNumber);
        m.def("ElectronegativityForElement", elec_for_element, "Allen-scale electronegativity");
        m.def("GetBondsAnglesDihedrals", get_bonds_angles_dihedrals);
        m.def("SelectionIsClosed", SelectionIsClosed, arg("mol"), arg("ids"), arg("structure_only")=false);

        class_<SmartsPattern>(m, "SmartsPattern")
            .def(init<std::string const&>())
            .def("atomCount", &SmartsPattern::atomCount)
            .def("pattern",   &SmartsPattern::pattern)
            .def("warnings",  &SmartsPattern::warnings)
            .def("findMatches", &SmartsPattern::findMatches)
            .def("match",     &SmartsPattern::match)
            ;
    }
}}
