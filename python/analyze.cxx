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
    void assign_1(SystemPtr mol, bool compute_resonant_charges, int timeout_ms) {
        unsigned flags = 0;
        if (compute_resonant_charges) flags |= AssignBondOrder::ComputeResonantCharges;
        AssignBondOrderAndFormalCharge(mol, flags, std::chrono::milliseconds(timeout_ms));
    }
    void assign_2(SystemPtr mol, list ids, bool compute_resonant_charges, int timeout_ms) {
        unsigned flags = 0;
        if (compute_resonant_charges) flags |= AssignBondOrder::ComputeResonantCharges;
        AssignBondOrderAndFormalCharge(mol, ids_from_python(ids), INT_MAX, flags, std::chrono::milliseconds(timeout_ms));
    }
    void assign_3(SystemPtr mol, list ids, int total_charge, bool compute_resonant_charges, int timeout_ms) {
        unsigned flags = 0;
        if (compute_resonant_charges) flags |= AssignBondOrder::ComputeResonantCharges;
        AssignBondOrderAndFormalCharge(mol, ids_from_python(ids), total_charge, flags, std::chrono::milliseconds(timeout_ms));
    }

    dict find_distinct_fragments(SystemPtr mol, object keys_obj) {
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        std::vector<std::string> keys;
        size_t i, n = len(keys_obj);
        for (i=0; i<n; i++) keys.push_back(extract<std::string>(keys_obj[i]));
        dict result;
        for (auto& iter : FindDistinctFragments(mol, fragments, keys)) {
            result[iter.first] = to_python(iter.second);
        }
        return result;
    }

    list find_matches(SmartsPattern const& s, AnnotatedSystem const& sys, 
                                              list const& starts) {
        return to_python(s.findMatches(sys, ids_from_python(starts)));
    }

    list wrap_sssr(SystemPtr mol, list atoms, bool all_relevant=false)
    {
        return to_python(GetSSSR(mol, ids_from_python(atoms), all_relevant));
    }
    list ring_systems(SystemPtr mol, list atoms) {
        return to_python(RingSystems(mol, GetSSSR(mol, ids_from_python(atoms), true)));
    }
    double elec_for_element(int n) {
        return DataForElement(n).eneg;
    }

    list compute_topids(SystemPtr mol) {
        return to_python(ComputeTopologicalIds(mol));
    }

    void guess_hydrogen(SystemPtr mol, list ids) {
        GuessHydrogenPositions(mol, ids_from_python(ids));
    }

    list convert_vector_idlist(std::vector<IdList> const& v) {
        list L;
        for (auto ids : v) {
            list s;
            for (auto id : ids) {
                s.append(id);
            }
            L.append(s);
        }
        return L;
    }

    dict get_bonds_angles_dihedrals(SystemPtr mol) {
        std::vector<IdList> non_pseudo_bonds;
        std::vector<IdList> pseudo_bonds;
        std::vector<IdList> angles;
        std::vector<IdList> dihedrals;
        GetBondsAnglesDihedrals(mol, mol->atoms(), non_pseudo_bonds, pseudo_bonds, angles, dihedrals);
        dict d;
        d["non_pseudo_bonds"] = convert_vector_idlist(non_pseudo_bonds);
        d["pseudo_bonds"] = convert_vector_idlist(pseudo_bonds);
        d["angles"] = convert_vector_idlist(angles);
        d["dihedrals"] = convert_vector_idlist(dihedrals);
        return d;
    }

}

namespace desres { namespace msys { 

    void export_analyze() {
        def("AssignBondOrderAndFormalCharge", assign_1);
        def("AssignBondOrderAndFormalCharge", assign_2);
        def("AssignBondOrderAndFormalCharge", assign_3);
        /* Yes, we have two interfaces for SSSR, this one and the one in
         * AnnotatedSystem.  This one lets you specify which atoms you
         * want the rings for, and doesn't force you to do any annotation, 
         * which is what we want.  AnnotatedSystem's rings() method only
         * lets you find rings connected to specific atoms or bonds. 
         */
        def("GetSSSR", wrap_sssr);
        def("RingSystems", ring_systems);
        def("ComputeTopologicalIds", compute_topids);
        def("GuessBondConnectivity", GuessBondConnectivity);
        def("FindDistinctFragments", find_distinct_fragments);
        def("RadiusForElement", RadiusForElement);
        def("MassForElement", MassForElement);
        def("PeriodForElement", PeriodForElement);
        def("GroupForElement", GroupForElement);
        def("ElementForAbbreviation", ElementForAbbreviation);
        def("GuessHydrogenPositions", guess_hydrogen);
        def("AbbreviationForElement", AbbreviationForElement);
        def("Analyze", Analyze);
        def("GuessAtomicNumber", GuessAtomicNumber);
        def("ElectronegativityForElement", elec_for_element, "Allen-scale electronegativity");
        def("GetBondsAnglesDihedrals", get_bonds_angles_dihedrals);

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

