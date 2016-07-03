#include "../analyze.hxx"
#include "../sssr.hxx"
#include "../elements.hxx"
#include <unordered_map>

namespace lpsolve {
#ifdef __APPLE__
#include <lp_lib.h>
#else
#include <lp_solve/lp_lib.h>
#endif
}

using namespace desres::msys;

struct range_t {
    int lb, ub;
    range_t() {}
    range_t(int&& l, int&& u) : lb(l), ub(u) {}
};

typedef std::unordered_map<Id,range_t> lpmap_t;

static int max_free_pairs(SystemPtr mol, const Id aid) {
    int maxFree;
    int nbonds=mol->bondCountForAtom(aid);
    int anum=mol->atomFAST(aid).atomic_number;
    int group=GroupForElement(anum);
    int period=PeriodForElement(anum);
    if (anum<3) {
        /* hydrogen and helium */
        maxFree=1-nbonds;
    } else if (group>=1 && group <=12) {
        /* metals / transition metals */
        maxFree=0;
    } else {
        /* Everything else */
        maxFree=4-nbonds;
        if (period>1 && nbonds>2) maxFree=std::max(maxFree+1, 1);
    }
    return std::max(0, maxFree);
}

static int maxbond_lp(SystemPtr mol, const Id aid0, const Id aid1){
    
    int maxOrders[2];
    Id ids[2]={aid0,aid1};
    for(int idx=0;idx<2;++idx){
        Id aid=ids[idx];  
        Id nbonds = mol->bondCountForAtom(aid);
        int anum=mol->atomFAST(aid).atomic_number;

        /* now determine max bond order */
        int maxbo;
        switch (anum) {
        case 1:  // H
        case 9:  // F
            // single bond max to hydrogen, flourine
            maxbo=1;
            break;
            // other halogens (per[chlor,brom,iod]ates)
        case 17: // Cl
        case 35: // Br
        case 53: // I
            if (nbonds==1) {
                maxbo=1;
            } else {
                maxbo=2;
            }
            break;
        case 5: // B
        case 6: // C
        case 7: // N
            if (nbonds == 4) {
                // single bond max to boron-, saturated carbon, nitrogen+
                maxbo=1;
            } else if(nbonds == 3) {
                // double bond max to boron-, carbon, nitrogen+ w/3 bonds
                maxbo=2;
            } else {
                maxbo=3;
            }
            break;
        case 8:  // O
            if (nbonds == 1){
                maxbo=3;
            } else {
                /* allows [O+](-R)(=R)  */
                maxbo=2;
            }
            break;

        default:
            /* Catch all... Should be fairly conservative */
            if (nbonds<3) {
                maxbo=3;
            } else if (nbonds<6) {
                maxbo=2;
            } else {
                maxbo=1;         
            }
        }
        maxOrders[idx]=maxbo;
    }
    return std::min(maxOrders[0],maxOrders[1]);
}

static bool allow_hextet_for_atom(SystemPtr mol, Id aid1) {
    int anum1=mol->atom(aid1).atomic_number;
    Id nbonds=mol->bondCountForAtom(aid1);
    /* Allow hextets for... */
    if ( 
        ( (anum1==5 || anum1==13) && (nbonds<4)               ) || // Al, B
        ( (anum1==6 || anum1==14) && (nbonds<4) ) || // unsaturated carbon
        ( (anum1==7 || anum1==15) && (nbonds<3) ) || // unsaturated nitrogen
        ( (anum1==8 || anum1==16) && (nbonds<2) )    // unsaturated oxygen
        ) return true;
    return false;
}
                
void desres::msys::AssignAromaticBondOrders(SystemPtr mol) {
    MultiIdList fragments;
    mol->updateFragids(&fragments);
    for (auto const& frag : fragments) {
        auto rings = GetSSSR(mol, frag, true);

        // consider only rings with <8 atoms and all atoms having <4 bonds
        MultiIdList keepRings;
        MultiIdList ringBonds;  // list of bond ids of rings
        for (auto& ring : rings) {
            if (ring.size()>7) continue;
            bool good = true;
            for (auto id : ring) {
                good = good && mol->bondCountForAtom(id)>3;
            }
            if (good) {
                IdList bonds;
                for (Id i=0; i<ring.size(); i++) {
                    bonds.push_back(mol->findBond(ring[i], ring[(i+1)%ring.size()]));
                }
                ringBonds.emplace_back(std::move(bonds));
                keepRings.emplace_back(std::move(ring));
            }
        }
        auto systems = RingSystems(mol, keepRings);
        for (auto& system : systems) {
            if (system.size()==1) continue;     // already have this one
            IdList bonds;
            for (Id rid : system) {
                for (Id bid : ringBonds[rid]) {
                    bonds.push_back(bid);
                }
            }
            sort_unique(bonds);
            ringBonds.emplace_back(std::move(bonds));
        }
        
        // build LPs
        lpmap_t atom_lp, bond_lp;

        int totalValence = 0;
        for (auto ai : frag) {
            int anum = mol->atomFAST(ai).atomic_number;
            auto const& data = DataForElement(anum);
            totalValence += data.nValence;
            atom_lp.emplace(ai, range_t(0,max_free_pairs(mol, anum)));
            for (auto b : mol->bondsForAtom(ai)) {
                auto aj = mol->bondFAST(b).other(ai);
                if (ai>aj) continue;
                bond_lp.emplace(b, range_t(1, maxbond_lp(mol, ai, aj)));
            }
        }

        /* presolve for unknown bond orders / lp counts.
         * can only presolve for cases where there is a single unknown.
         */
        
        auto unsolved=frag;
        IdList keep;
        keep.reserve(frag.size());
    
        int npresolve=0;
        int nsolved=0;
    
        bool stillsolving;
        do {
            stillsolving=false;
            for (auto aid1 : unsolved) {
                atom_t const& atm1=mol->atomFAST(aid1);
                int period=PeriodForElement(atm1.atomic_number);
                if(period>2){
                    keep.push_back(aid1);
                    continue;
                }
                int octval=period<=1 ? 1 : 4; 
                
                int unkcount=0;
                auto lastatom=atom_lp.end();
                auto lastbond=bond_lp.end();
                
                auto aiter=atom_lp.find(aid1);
                if(aiter->second.lb==aiter->second.ub){
                    octval-=aiter->second.lb;
                }else{
                    unkcount++;
                    lastatom=aiter;
                } 
                
                auto bonds = mol->bondsForAtom(aid1);
                for (auto bid : bonds) {
                    auto biter=bond_lp.find(bid);
                    if (biter->second.lb==biter->second.ub){
                        octval-=biter->second.lb;
                    }else{
                        unkcount++;
                        lastbond=biter;
                    }
                }
                if(unkcount==1){
                    // Get aid/bid that we can now solve for
                    if( lastatom != atom_lp.end()){
                        auto &fixrange= lastatom->second;
    
                        if(!(octval >=fixrange.lb && octval <= fixrange.ub) && bonds.size()>0) {
                            MSYS_FAIL("Valence count solution is invalid for atom "<< lastatom->first <<
                                " ( " << octval << " " << fixrange.ub <<" )");
                        }
    
                        /* Only solve if possible lp count==0 or cant form hextet */
                        if(octval==0 || !allow_hextet_for_atom(mol, aid1) ){
                            fixrange.lb=octval;
                            fixrange.ub=octval;
                            
                            stillsolving=true;
                            nsolved+=1;
                        } else {
                            keep.push_back(aid1);
                        }
                    } else {
                        auto &fixrange= lastbond->second;
                        if(!(octval >=fixrange.lb && octval <= fixrange.ub)) {
                            MSYS_FAIL("Bond order solution is invalid");
                        }
                        fixrange.lb=octval;
                        fixrange.ub=octval;
    
                        stillsolving=true;
                        nsolved+=1;
                    }
                }else if (unkcount > 1){
                    keep.push_back(aid1);
                }
            }
            unsolved.clear();
            unsolved.swap(keep);
            npresolve++;
        } while (stillsolving);
    }
}

