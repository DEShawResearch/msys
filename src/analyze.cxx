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

    int FragmentChargeFromAtomicCharge(SystemPtr mol,
                                        IdList const& atoms){
        
        double qSum=0;
        double qMax=0;
        BOOST_FOREACH(Id a, atoms) {
            double aq=mol->atom(a).charge;
            qSum+=aq;
            aq=fabs(aq);
            if(aq>qMax) qMax=aq;
        }

        if(qMax!=0) return int(qSum);
        return INT_MAX;
    }


    int FragmentChargeFromFormalCharge(SystemPtr mol,
                                       IdList const& atoms){
        
        int fcSum=0;
        int fcMax=0;
        BOOST_FOREACH(Id a, atoms) {
            int fq=mol->atom(a).formal_charge;
            fcSum+=fq;
            fq=abs(fq);
            if(fq>fcMax) fcMax=fq;
        }

        if(fcMax) return fcSum;
        return INT_MAX;
    }

    int FragmentChargeFromBondOrders(SystemPtr mol,
                                     IdList const& atoms){
        
        int bqSum=0;
        int bqMax=0;
        bool hasBO=false;
        bool fullOctet=true;
        BOOST_FOREACH(Id a, atoms) {
            atom_t const& atom = mol->atom(a);
            int anum=atom.atomic_number;
            int group=GroupForElement(anum);

            IdList bonds = mol->bondsForAtom(a);

            int boSum=0;
            bool hasPolar=false;
            bool hasN4=false;
            if(anum!=6){
                BOOST_FOREACH(Id b, bonds) {
                    bond_t const& bond = mol->bond(b);
                    boSum+=bond.order;
                }
            }else{
                std::vector<Id> otherN;
                BOOST_FOREACH(Id b, bonds) {
                    bond_t const& bond = mol->bond(b);
                    boSum+=bond.order;
                    Id a2=bond.other(a);
                    int anum2=mol->atom(a2).atomic_number;
                    if(anum2==7){
                        hasPolar=true;
                        otherN.push_back(a2);
                    }else if(anum2==8 || anum2==16 || anum2==9 || anum2==17){
                        hasPolar=true;
                    }
                }
                if(boSum==3){
                    BOOST_FOREACH(Id a2, otherN){
                        int boSum2=0;
                        BOOST_FOREACH(Id b, mol->bondsForAtom(a2)) {
                            boSum2+= mol->bond(b).order;
                        }
                        if(boSum2==4){
                            hasN4=true;
                            break;
                        }
                    }
                }
            }

            int nbonds=bonds.size();
            hasBO |= (boSum!=nbonds);
            
            int q=0;
            if(anum==1 && boSum!=1){
                q=1;
            }else if( (anum==5 || anum==13) && boSum==4){
                q=1;
            }else if(anum==6 && boSum==3){
                if(nbonds==1 || hasN4){
                    q=-1;
                }else{
                    q= hasPolar ? 1 : -1;
                }
            }else if(anum==7){
                if(boSum==2){
                    q=-1;
                }else if(boSum==4){
                    q=1;
                }
            }else if(anum==8){
                if(boSum==1){
                    q=-1;
                }else if(boSum==3){
                    q=1;
                }
            }else if(anum==15){
                if(boSum==4){
                    q=1;
                }else if(boSum==2){
                    q=-1;
                }
            }else if(anum==16){
                if(boSum==1 || boSum==5){
                    q=-1;
                }else if(boSum==3){
                    q=1;
                }
            }else if(anum==17){
                if(boSum==0){
                    q=-1;
                }else if(boSum==4){
                    q=3;
                }
            }else if(boSum==0){
                if(group==1){
                    q=1;
                }else if(group==2 || anum==30){
                    q=2;
                }else if(group==17){
                    q=-1;
                }
            }

            int electrons = DataForElement(anum).nValence - boSum - q;
            if (electrons < 0 || electrons % 2){
                fullOctet=false;
            }else if(anum==1){
                fullOctet &= (2*boSum+electrons==2);
            }else{
                fullOctet &= (2*boSum+electrons==8);               
            }

            bqSum+=q;
            q=abs(q);
            if(q>bqMax)bqMax=q;
        }

        if(hasBO || fullOctet) return bqSum;
        return INT_MAX;
    }


    int GuessBestFragmentCharge(SystemPtr mol,
                                IdList const& atoms){

        int qTot=FragmentChargeFromFormalCharge(mol, atoms);
        if(qTot!=INT_MAX) return qTot;
        qTot=FragmentChargeFromAtomicCharge(mol, atoms);
        if(qTot!=INT_MAX) return qTot;        
        return FragmentChargeFromBondOrders(mol, atoms);
    }

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
                            def_bond, def_angle, -M_PI/2+def_dihedral);
                    hyd1.x = pos.x;
                    hyd1.y = pos.y;
                    hyd1.z = pos.z;

                    pos = apply_dihedral_geometry(
                            Vec3(&mol->atom(A).x),
                            Vec3(&mol->atom(C).x),
                            Vec3(&mol->atom(r).x),
                            def_bond, def_angle, -M_PI/2-def_dihedral);
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
