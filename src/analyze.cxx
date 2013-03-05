#include "analyze.hxx"
#include "analyze/bond_orders.hxx"
#include "elements.hxx"
#include "graph.hxx"
#include "geom.hxx"
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

    IdList FindDistinctFragments(SystemPtr mol, MultiIdList const& fragments) {
        IdList result;
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

    static const double def_bond = 1.0;
    static const double def_angle = 109.5 * M_PI/180;
    static const double def_dihedral = -120.0 * M_PI/180;

    void GuessHydrogenPositions(SystemPtr mol, IdList const& hatoms) {
        /* partition by root atom */
        typedef std::map<Id,IdList> RootMap;
        RootMap map;
        BOOST_FOREACH(Id h, hatoms) {
            if (mol->bondCountForAtom(h)!=1) continue;
            map[mol->bondedAtoms(h).at(0)].push_back(h);
        }

        /* process each root atom */
        for (RootMap::iterator it=map.begin(); it!=map.end(); ++it) {
            const Id r = it->first;
            const atom_t& root = mol->atom(r);
            IdList& hlist = it->second;
            std::sort(hlist.begin(), hlist.end());

            /* get other atoms bonded to root. */
            IdList c;
            BOOST_FOREACH(Id b, mol->bondedAtoms(r)) {
                if (!std::binary_search(hlist.begin(), hlist.end(), b)) {
                    c.push_back(b);
                }
            }

            if (hlist.size()==2) {
                atom_t& hyd1 = mol->atom(hlist[0]);
                atom_t& hyd2 = mol->atom(hlist[1]);
                if (c.empty()) {
                    /* water */
                    hyd1.x = root.x;
                    hyd1.y = root.y;
                    hyd1.z = root.z + def_bond;
                    Vec3 pos = apply_dihedral_geometry(
                            Vec3(&mol->atom(hlist[0]).x),
                            Vec3(&mol->atom(hlist[0]).x),
                            Vec3(&mol->atom(r).x),
                            def_bond, def_angle, 0);
                    hyd2.x = pos.x;
                    hyd2.y = pos.y;
                    hyd2.z = pos.z;

                } else if (c.size()==1) {
                    Id C = c[0];
                    Id A = C;
                    if (mol->bondCountForAtom(C)>1) {
                        BOOST_FOREACH(A, mol->bondedAtoms(C)) {
                            if (A!=r) break;
                        }
                    }
                    Vec3 pos = apply_dihedral_geometry(
                            Vec3(&mol->atom(A).x),
                            Vec3(&mol->atom(C).x),
                            Vec3(&mol->atom(r).x),
                            def_bond, def_angle, def_dihedral);
                    hyd1.x = pos.x;
                    hyd1.y = pos.y;
                    hyd1.z = pos.z;

                    pos = apply_dihedral_geometry(
                            Vec3(&mol->atom(A).x),
                            Vec3(&mol->atom(C).x),
                            Vec3(&mol->atom(r).x),
                            def_bond, def_angle, -def_dihedral);
                    hyd2.x = pos.x;
                    hyd2.y = pos.y;
                    hyd2.z = pos.z;
                }

            } else if (hlist.size()==1) {
                atom_t& hyd1 = mol->atom(hlist[0]);

                if (c.empty()) {
                    hyd1.x = root.x;
                    hyd1.y = root.y;
                    hyd1.z = root.z + def_bond;

                } else if (c.size()==1) {
                    /* find another atom bonded to c to define the plane */
                    Id C = c[0];
                    Id A = C;
                    if (mol->bondCountForAtom(C)>1) {
                        BOOST_FOREACH(A, mol->bondedAtoms(C)) {
                            if (A!=r) break;
                        }
                    }
                    Vec3 pos = apply_dihedral_geometry(
                            Vec3(&mol->atom(A).x),
                            Vec3(&mol->atom(C).x),
                            Vec3(&mol->atom(r).x),
                            def_bond, def_angle, 0);
                    hyd1.x = pos.x;
                    hyd1.y = pos.y;
                    hyd1.z = pos.z;

                } else {
                    Float ux=0, uy=0, uz=0;
                    BOOST_FOREACH(Id b, c) {
                        ux += mol->atom(r).x - mol->atom(b).x;
                        uy += mol->atom(r).y - mol->atom(b).y;
                        uz += mol->atom(r).z - mol->atom(b).z;
                    }
                    Float len = sqrt(ux*ux + uy*uy + uz*uz);
                    if (len) {
                        ux /= len;
                        uy /= len;
                        uz /= len;
                    }
                    hyd1.x = root.x + def_bond*ux;
                    hyd1.y = root.y + def_bond*uy;
                    hyd1.z = root.z + def_bond*uz;
                }
            }
        }
    }
}}
