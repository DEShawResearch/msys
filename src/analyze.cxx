#include "analyze.hxx"
#include "analyze/bond_orders.hxx"

namespace desres { namespace msys {

    void AssignBondOrderAndFormalCharge(SystemPtr mol) {
        MultiIdList frags;
        mol->updateFragids(&frags);
        BOOST_FOREACH(IdList const& frag, frags) {
            AssignBondOrderAndFormalCharge(mol, frag);
        }
    }
    
    void AssignBondOrderAndFormalCharge(SystemPtr mol, 
                                        IdList const& atoms,
                                        int total_charge) {
#ifdef MSYS_WITHOUT_LPSOLVE
        MSYS_FAIL("LPSOLVE functionality was not included.");
#else
        if (atoms.empty()) return;
        BondOrderAssignerPtr boa=BondOrderAssigner::create(mol, atoms);
        if (total_charge != INT_MAX) {
            boa->setTotalCharge(total_charge);
        }
        boa->solveIntegerLinearProgram();
        boa->assignSolutionToAtoms();
#endif
    }

}}
