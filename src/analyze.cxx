#include "analyze.hxx"
#include "analyze/bond_orders.hxx"
#include "elements.hxx"
#include <periodicfix/contacts.hxx>
#include <stdio.h>

namespace {
    using namespace desres::msys;
    struct BondFinder {
        SystemPtr mol;
        BondFinder(SystemPtr m) : mol(m) {}

        bool exclude(Id i, Id j) {
            if (i>=j) return true;
            int ai = mol->atom(i).atomic_number;
            int aj = mol->atom(j).atomic_number;
            if ((ai==1 && aj==1) ||
                (ai==0 && aj==0)) 
                return true;
            return false;
        }

        void operator()(Id i, Id j, double d2) {
            if (d2>0.001) {
                int ai = mol->atom(i).atomic_number;
                int aj = mol->atom(j).atomic_number;
                double ri = RadiusForElement(ai);
                double rj = RadiusForElement(aj);
                double cut = 0.6 * (ri+rj);
                if (d2 < cut*cut) {
                    mol->addBond(i,j);
                }
            }
        }
    };
}

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


    void GuessBondConnectivity(SystemPtr mol) {
        std::vector<Float> pos(3*mol->maxAtomId());
        if (pos.empty()) return;
        IdList atoms(mol->atoms());
        BOOST_FOREACH(Id i, atoms) {
            atom_t const& atom = mol->atom(i);
            pos[3*i  ] = atom.x;
            pos[3*i+1] = atom.y;
            pos[3*i+2] = atom.z;
        }
        BondFinder finder(mol);
        periodicfix::find_contacts(4.0, &pos[0],
                                   atoms.begin(), atoms.end(),
                                   atoms.begin(), atoms.end(),
                                   finder);
    }
}}
