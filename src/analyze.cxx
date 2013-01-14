#include "analyze.hxx"
#include "analyze/bond_orders.hxx"
#include "elements.hxx"
#include "graph.hxx"
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
        /* if a hydrogen has multiple bonds, keep the shortest one */
        for (Id i=0; i<mol->maxAtomId(); i++) {
            if (!mol->hasAtom(i)) continue;
            if (mol->atom(i).atomic_number!=1) continue;
            if (mol->bondCountForAtom(i)<=1) continue;
            Id shortest_bond = BadId;
            double shortest_dist = HUGE_VAL;
            const double* pi = &pos[3*i];
            const double x=pi[0];
            const double y=pi[1];
            const double z=pi[2];
            IdList bonds = mol->bondsForAtom(i);
            BOOST_FOREACH(Id b, bonds) {
                Id j = mol->bond(b).other(i);
                if (j>i) continue;
                const double* pj = &pos[3*j];
                const double dx = pj[0]-x;
                const double dy = pj[1]-y;
                const double dz = pj[2]-z;
                const double d2 = dx*dx + dy*dy + dz*dz;
                if (d2<shortest_dist) {
                    shortest_bond = b;
                    shortest_dist = d2;
                }
            }
            assert(!bad(shortest_bond));
            BOOST_FOREACH(Id b, bonds) {
                if (b!=shortest_bond) mol->delBond(b);
            }
        }
    }

    IdList FindDistinctFragments(SystemPtr mol) {
        IdList result;
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        /* will compute graphs lazily */
        std::vector<GraphPtr> graphs(fragments.size());
        typedef std::map<std::string, IdList> FragmentHash;
        FragmentHash fragment_hash;
        for (Id i=0; i<fragments.size(); i++) {
            fragment_hash[Graph::hash(mol, fragments[i])].push_back(i);
        }
        FragmentHash::iterator it;
        for (it=fragment_hash.begin(); it!=fragment_hash.end(); ++it) {
            /* unique formula -> unique fragment */
            IdList& frags = it->second;
            if (frags.size()==1) {
                result.push_back(frags[0]);
                continue;
            }
            /* must do isomorphism checks. */
            BOOST_FOREACH(Id frag, frags) {
                graphs[frag] = Graph::create(mol, fragments[frag]);
            }
            std::vector<IdPair> perm;
            while (!frags.empty()) {
                result.push_back(frags[0]);
                IdList unmatched;
                GraphPtr ref = graphs[frags[0]];
                for (Id i=1; i<frags.size(); i++) {
                    GraphPtr sel = graphs[frags[i]];
                    if (!ref->match(sel, perm)) {
                        unmatched.push_back(frags[i]);
                    }
                }
                frags.swap(unmatched);
            }
        }
        std::sort(result.begin(), result.end());
        return result;
    }
}}
