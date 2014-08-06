#include "analyze.hxx"
#include "analyze/bond_orders.hxx"
#include "elements.hxx"
#include "graph.hxx"
#include "geom.hxx"
#include "contacts.hxx"
#include <stdio.h>
#include <boost/foreach.hpp>

namespace {
    using namespace desres::msys;
    struct BondFinder {
        SystemPtr mol;
        BondFinder(SystemPtr m) : mol(m) {}

        bool exclude(Id i, Id j) const {
            if (i>=j) return true;
            int ai = mol->atomFAST(i).atomic_number;
            int aj = mol->atomFAST(j).atomic_number;
            if ((ai==1 && aj==1) ||
                (ai==0 && aj==0)) 
                return true;
            return false;
        }

        void operator()(Id i, Id j, double d2) const {
            if (d2>0.001) {
                int ai = mol->atomFAST(i).atomic_number;
                int aj = mol->atomFAST(j).atomic_number;
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
        MultiIdList fragments;
        mol->updateFragids(&fragments);
        IdList pmap(mol->maxAtomId(), BadId);

        /* will compute graphs lazily */
        std::vector<GraphPtr> graphs(fragments.size());
        typedef std::map<std::string, IdList> FragmentHash;
        FragmentHash fragment_hash;

        for (Id i=0; i<fragments.size(); i++) {
            fragment_hash[Graph::hash(mol, fragments[i])].push_back(i);
        }
        FragmentHash::iterator it;
        for (it=fragment_hash.begin(); it!=fragment_hash.end(); ++it) {
            /* fetch fragments with the same formula */
            IdList& frags = it->second;
            /* unique formula -> unique fragment */
            if (frags.size()==1) {
                AssignBondOrderAndFormalCharge(mol, fragments[frags[0]]);
                continue;
            }

            /* We have multiple fragments with the same formula.  */
            BOOST_FOREACH(Id frag, frags) {
                graphs[frag] = Graph::create(mol, fragments[frag]);
            }
            std::vector<IdPair> perm;
            while (!frags.empty()) {
                AssignBondOrderAndFormalCharge(mol, fragments[frags[0]]);
                IdList unmatched;
                GraphPtr ref = graphs[frags[0]];
                for (Id i=1; i<frags.size(); i++) {
                    GraphPtr sel = graphs[frags[i]];
                    if (!ref->match(sel, perm)) {
                        /* didn't match, so push into the next iteration 
                         * for another bond order calculation */
                        unmatched.push_back(frags[i]);
                    } else {
                        /* map atom properties */
                        BOOST_FOREACH(IdPair const&p, perm) {
                            const Id ai = p.first;
                            const Id bi = p.second;
                            mol->atom(bi).formal_charge = mol->atom(ai).formal_charge;
                            mol->atom(bi).resonant_charge = mol->atom(ai).resonant_charge;
                            pmap.at(ai) = bi;
                        }
                        /* map bond properties */
                        BOOST_FOREACH(IdPair const&p, perm) {
                            const Id ai = p.first;
                            const Id bi = p.second;
                            BOOST_FOREACH(Id bnd, mol->bondsForAtom(ai)) {
                                bond_t const& src = mol->bond(bnd);
                                const Id aj = src.other(ai);
                                if (ai>aj) continue;
                                const Id bj = pmap.at(aj);
                                bond_t& dst = mol->bond(mol->findBond(bi,bj));
                                dst.order = src.order;
                                dst.resonant_order = src.resonant_order;
                            }
                        }
                    }
                }
                frags.swap(unmatched);
            }
        }
    }

#ifdef MSYS_WITHOUT_LPSOLVE
    IdList ComputeTopologicalIds(SystemPtr mol) {
        MSYS_FAIL("LPSOLVE functionality was not included.");
        return IdList();
    }
#else
    // defined in src/analyze/topological_ids.cxx
#endif
    
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
        static const double cutoff = 4.0;
        if (pos.empty()) return;
        IdList atoms(mol->atoms());
        BOOST_FOREACH(Id i, atoms) {
            atom_t const& atom = mol->atom(i);
            pos[3*i  ] = atom.x;
            pos[3*i+1] = atom.y;
            pos[3*i+2] = atom.z;
        }
        BondFinder finder(mol);
        if (mol->ctCount()==1) {
            find_contacts(cutoff, &pos[0],
                          atoms.begin(), atoms.end(),
                          atoms.begin(), atoms.end(),
                          finder);
        } else {
            for (Id i=0, n=mol->ctCount(); i<n; i++) {
                IdList const& atoms = mol->atomsForCt(i);
                find_contacts(cutoff, &pos[0],
                              atoms.begin(), atoms.end(),
                              atoms.begin(), atoms.end(),
                              finder);
            }
        }
        /* if a hydrogen has multiple bonds, keep the shortest one */
        for (Id i=0; i<mol->maxAtomId(); i++) {
            if (!mol->hasAtom(i)) continue;
            if (mol->atomFAST(i).atomic_number!=1) continue;
            if (mol->bondCountForAtom(i)<=1) continue;
            Id shortest_bond = BadId;
            double shortest_dist = HUGE_VAL;
            const double* pi = &pos[3*i];
            const double x=pi[0];
            const double y=pi[1];
            const double z=pi[2];
            IdList bonds = mol->bondsForAtom(i);    /* yes, make a copy! */
            BOOST_FOREACH(Id b, bonds) {
                Id j = mol->bond(b).other(i);
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

    IdList AddHydrogens(SystemPtr mol, IdList const &parents) {
        IdList added;

        BOOST_FOREACH(Id atm, parents) {
            int anum = mol->atom(atm).atomic_number;
            int fq = mol->atom(atm).formal_charge;
            int grp = GroupForElement(anum);
            if (anum <= 2 || grp < 13) continue;

            int degree = 0;
            int valence = 0;
            BOOST_FOREACH(Id bnd, mol->bondsForAtom(atm)) {
                ++degree;
                valence += mol->bond(bnd).order;
            }
            /* shortcut check for saturated atoms */
            if( (degree>3 || valence>3) && PeriodForElement(anum)<3){
                continue;
            }

            int target = grp-10;
            int electrons = target - valence - fq;
            int extraH = 0;
            
            // add protons to get rid of radical
            if (electrons % 2) {
                ++extraH;
                --electrons;
                ++degree;
                ++valence;
            }

            int nlp = electrons/2; // number of lone pairs
            /* convert lone pairs to 2x hydrogens to satisfy octets.
             * This may not hold for hypervalent species
             */
            int more = 4 - nlp - valence;
            extraH += std::max(0, 2*std::min(more, nlp));

            for (int i=0; i<extraH; i++) {
                Id H = mol->addAtom(mol->atom(atm).residue);
                mol->addBond(atm,H);
                mol->atom(H).fragid = mol->atom(atm).fragid;
                mol->atom(H).atomic_number = 1;
                mol->atom(H).name = "H";
                added.push_back(H);
            }
        }
        return added;
    }

    static const double def_bond = 1.0;
    static const double def_angle = 109.5 * M_PI/180;
    static const double def_dihedral = -120.0 * M_PI/180;

    void GuessHydrogenPositions(SystemPtr mol, IdList const& hatoms) {

        /* Randomize orientation of water molecules, deterministically. */
        srand48(1999);

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

                    /* sample a point on the unit sphere */
                    double u1 = drand48();
                    double u2 = drand48();
                    double z = 2*u1 - 1.0;
                    double phi = 2*M_PI*u2;
                    double R = sqrt(1-z*z);

                    /* orient and place the first hydrogen */
                    hyd1.x = root.x + def_bond * R * cos(phi);
                    hyd1.y = root.y + def_bond * R * sin(phi);
                    hyd1.z = root.z + def_bond * z;

                    /* arbitrary reference position for the second hydrogen */
                    Float A[3] = {hyd1.x, hyd1.y + def_bond, hyd1.z};
                    /* position with random rotation about O-H1 */
                    apply_dihedral_geometry(hyd2.pos(),
                            A,
                            mol->atom(hlist[0]).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, 2*M_PI*drand48());

                } else if (c.size()==1) {
                    Id C = c[0];
                    Id A = C;
                    if (mol->bondCountForAtom(C)>1) {
                        BOOST_FOREACH(A, mol->bondedAtoms(C)) {
                            if (A!=r) break;
                        }
                    }
                    apply_dihedral_geometry(hyd1.pos(),
                            mol->atom(A).pos(),
                            mol->atom(C).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, -M_PI/2+def_dihedral);

                    apply_dihedral_geometry(hyd2.pos(),
                            mol->atom(A).pos(),
                            mol->atom(C).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, -M_PI/2-def_dihedral);
                } else if (c.size()==2) {
                    apply_dihedral_geometry(hyd1.pos(),
                            mol->atom(c[0]).pos(),
                            mol->atom(c[1]).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, def_dihedral);
                    apply_dihedral_geometry(hyd2.pos(),
                            mol->atom(c[1]).pos(),
                            mol->atom(c[0]).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, def_dihedral);
                }

            } else if (hlist.size()==3) {
                if (c.size()>0) {
                    Id C = c[0];
                    Id A = C;
                    if (mol->bondCountForAtom(C)>1) {
                        BOOST_FOREACH(A, mol->bondedAtoms(C)) {
                            if (A!=r) break;
                        }
                    }
                    apply_dihedral_geometry(mol->atom(hlist[0]).pos(),
                            mol->atom(A).pos(),
                            mol->atom(C).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, def_dihedral/2);
                    apply_dihedral_geometry(mol->atom(hlist[1]).pos(),
                            mol->atom(hlist[0]).pos(),
                            mol->atom(C).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, def_dihedral);
                    apply_dihedral_geometry(mol->atom(hlist[2]).pos(),
                            mol->atom(C).pos(),
                            mol->atom(hlist[0]).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, def_dihedral);

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
                    apply_dihedral_geometry(hyd1.pos(),
                            mol->atom(A).pos(),
                            mol->atom(C).pos(),
                            mol->atom(r).pos(),
                            def_bond, def_angle, 0);

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
