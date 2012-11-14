#include <algorithm>  // sort
#include <cmath>     // fabs, pow, round
#include <stdexcept>
#include <sstream>
#include <cstdio>
#include <deque>
#include <boost/foreach.hpp>

#include "bond_orders.hxx"
#include "../sssr.hxx"
#include "aromatic.hxx"
#include "get_fragments.hxx"
#include "../quadsum.hxx"
#include "eigensystem.hxx"
#include "../elements.hxx"
#include "filtered_bonds.hxx"

using namespace desres::msys;

namespace lpsolve {
    /* for solving the integer linear equations. We put this in a new 
       namespace since it dosent already use one */ 
#include "lp_solve/lp_lib.h"
}


namespace {
    /* simple helper function to extract a value from a map (and assert that its found) */
    template<typename C>
    typename C::mapped_type& asserted_find(C &container, typename C::key_type const& key){
        typename C::iterator iter=container.find(key);
        if(iter == container.end()){
            throw std::runtime_error("Key not found in container");
        }
        return iter->second;
    }
   
    template<typename C>
    typename C::mapped_type const& asserted_find(C const& container, typename C::key_type const& key){
        typename C::const_iterator iter=container.find(key);
        if(iter == container.end()){
            throw std::runtime_error("Key not found in container");
        }
        return iter->second;
    }

    template<typename C>
    void print_Map(C const& m, std::string const& s){
        BOOST_FOREACH(typename C::value_type const& p,m){
            printf("%s: %u %d\n",s.c_str(),p.first,p.second); 
        }
    }
}


namespace desres { namespace msys {

    BondOrderAssignerPtr BondOrderAssigner::create(SystemPtr sys,  
                                                   IdList const& fragment){

        BondOrderAssignerPtr boa(new BondOrderAssigner);
        boa->_needRebuild=true;
        boa->_valid=false;
        boa->_total_charge_set=false;
        boa->_total_charge=0;
        boa->_mol=sys;

        BOOST_FOREACH(Id aid, fragment){
            assert(boa->_mol->hasAtom(aid));
            boa->_mol->atom(aid).formal_charge=0;
            boa->_mol->atom(aid).resonant_charge=0.0;
            BOOST_FOREACH(Id bid, boa->_mol->bondsForAtom(aid)){
                boa->_mol->bond(bid).order=1;
                boa->_mol->bond(bid).resonant_order=1.0;
            }
            if(boa->_mol->atom(aid).atomic_number<1){
                printf("BondOrderAssigner::create - Skipping atomid=%u with atomic_number<1\n",aid);
                continue;
            }
            boa->_fragatoms.push_back(aid);
        }

        MultiIdList allRings = GetSSSR(
                boa->_mol, boa->_fragatoms, true);
        /* FIXME: Add generator for fused rings <10 atoms from the above ssr ring set */ 
        
        BOOST_FOREACH(IdList & ring, allRings){
            std::set<Id> planar_atoms;
            BOOST_FOREACH(Id  aid1, ring){
                IdList bonded=filteredBondedAtoms(boa->_mol, aid1);
                planar_atoms.insert(bonded.begin(),bonded.end());
            }
            std::vector<Id> check(planar_atoms.begin(),planar_atoms.end());
            double planarity=ComputeRingPlanarity(boa->_mol, check);
            /* Algorithms in the bond order assigner are simplified if ring cycle 
               is closed, ie ring[0]==ring[-1] */
            size_t nratoms=ring.size();
            if(nratoms>1 && ring[0]!=ring[nratoms-1]) ring.push_back(ring[0]);
            if(planarity > boa->planarity_tolerance){
                boa->_nonplanar_rings.push_back(RingPair());
                RingPair &last=boa->_nonplanar_rings.back();
                last.first=planarity;
                swap(last.second,ring);
            }else{
                boa->_planar_rings.push_back(RingPair());
                RingPair &last=boa->_planar_rings.back();
                last.first=planarity;
                swap(last.second,ring);
            }
        }

        boa->_totalValence=0;
        electronRange tmprange;
        BOOST_FOREACH(Id aid1,boa->_fragatoms){
            atom_t& atm1=boa->_mol->atom(aid1);
            int anum1=atm1.atomic_number;
#if 0
            if(anum1>maxAtomicNum){
                std::stringstream msg;
                msg << "Atomic Number " << anum1 << " exceeds maximum allowed (" << maxAtomicNum << ")";
                throw std::runtime_error(msg.str());
            }

            IdList bonds=filteredBondsForAtom(boa->_mol,aid1);
            ChemData const& adata = DataForElement(anum1);
            if(adata==NODATA){
                printf("Warning: No property information available for Atomic Number %d",anum1);
                /* make sure missing data is for an ion */
                assert(bonds.size()==0);
            }
#else
            ChemData const& adata = DataForElement(anum1);
            IdList bonds=filteredBondsForAtom(boa->_mol,aid1);
#endif

            boa->_totalValence+=adata.nValence;

            tmprange.lb=0;
            tmprange.ub=adata.maxFree;
            boa->_atom_lp.insert(std::make_pair(aid1,tmprange));

            BOOST_FOREACH(Id bid, bonds){
                Id aid2 = boa->_mol->bond(bid).other(aid1);
                if(aid1>aid2) continue;

                int maxbo1=boa->max_bond_order(aid1);
                int maxbo2=boa->max_bond_order(aid2);
                if(maxbo1>maxbo2) maxbo1=maxbo2;

                tmprange.lb=1;
                tmprange.ub=maxbo1;
                boa->_bond_order.insert(std::make_pair(bid,tmprange));
            }
        }

        IdList unsolved;
        boa->presolve_octets(unsolved);
        if(unsolved.size()!=0){
            MultiIdList components;
            get_fragments(boa->_mol, unsolved, components);
#if DEBUGPRINT
            printf("Found %zu seperable components in fragment\n",components.size());
#endif
            for(Id cid=0;cid<components.size();++cid){
                boa->_component_assigners.push_back(ComponentAssigner::create(boa, components[cid], cid));
            }
        }

        /* Fill in bondinfo with presolved values */
        BOOST_FOREACH(electronMap::value_type const& epair, boa->_bond_order){
            int val=epair.second.lb;
            if(val!=epair.second.ub) continue;
            boa->_bondinfo.insert(solutionMap::value_type(epair.first, solutionValues(val,val)));
        }

        /* Fill in atominfo/chargeinfo with presolved values */
        boa->_presolved_charge=0;
        BOOST_FOREACH(electronMap::value_type const& epair, boa->_atom_lp){
            int val=epair.second.lb;
            if(val!=epair.second.ub) continue;
            Id aid=epair.first;
            boa->_atominfo.insert(solutionMap::value_type(aid, solutionValues(val,val)));
            /* Formal Charges are given by: fc[i]= ValenceElectrons[i] - freeElectrons[i] - 0.5*Sum_j ( BondElectrons[j] )
               'val' is free electron PAIRS, so multiply by 2 */
            int qtot=DataForElement(boa->_mol->atom(aid).atomic_number).nValence - 2*val;
            BOOST_FOREACH(Id const& bid, filteredBondsForAtom(boa->_mol,aid)){
                /* This should always be possible given how we determine atom lp data.
                 * If we know the lp's we also know all the bonds
                 * 'nonresonant' is the bond electron PAIRS, so no factor of 0.5 needed */
                qtot-=asserted_find(boa->_bondinfo, bid).nonresonant;
            }
            boa->_chargeinfo.insert(solutionMap::value_type(epair.first, solutionValues(qtot,qtot)));
            boa->_presolved_charge+=qtot;
        }
        
#if DEBUGPRINT
        printf("AfterPresolve: total_valence=%d for %zu atoms (%zu finalized,  q_of_finalized=%d)\n",
               boa->_totalValence, boa->_fragatoms.size(), boa->_atominfo.size(), boa->_presolved_charge);

#endif

        return boa;
   }

    void BondOrderAssigner::rebuild(){
        RingList new_nonplanar;
        RingList new_planar;
        BOOST_FOREACH(RingPair &rp, _nonplanar_rings){
            if(rp.first > planarity_tolerance){
                new_nonplanar.push_back(rp);
            }else{
                new_planar.push_back(rp);
            }
        }
        BOOST_FOREACH(RingPair &rp, _planar_rings){
            if(rp.first > planarity_tolerance){
                new_nonplanar.push_back(rp);
            }else{
                new_planar.push_back(rp);
            }
        }
        _nonplanar_rings.swap(new_nonplanar);
        _planar_rings.swap(new_planar);

        BOOST_FOREACH(ComponentAssignerPtr ca, _component_assigners){
            ca->build_integer_linear_program();
        }

        _needRebuild=false;
    }


    /* This seemingly simple function is important for 2 things... */
    int BondOrderAssigner::max_bond_order(const Id aid){
        /* default maximum bond order allowed per bond */
        static const int defmax=3;
        int maxbo;
        Id nbonds=filteredBondCountForAtom(_mol,aid);
        switch (_mol->atom(aid).atomic_number){
            case 1:  // H
            case 9:  // F
            case 17: // Cl
            case 35: // Br
            case 53: // I
                // single bond max to hydrogen, halogens
                maxbo=1;
                break;
            case 6: // C
            case 7: // N
                if(nbonds == 4){
                    // single bond max to saturated carbon, Nitrogen+
                    maxbo=1;
                }else if(nbonds == 3){
                    // double bond max to carbon, Nitrogen+ w/3 bonds
                    maxbo=2;
                }else{
                    maxbo=defmax;
                }
                break;
            case 8: // O
                /* double bond max for all oxygens, allows [O+](-R)(=R)  */
                maxbo=2;
                break;
            default:
                maxbo=defmax;
        }
        return maxbo;
    }


    /* Function that trys to presolve for unknown bond orders / lp counts 
       Its can only presolve for cases where there is a single unknown.
    */
    void BondOrderAssigner::presolve_octets(IdList &unsolved){
        //static desres::profiler::Symbol _("BondOrderAssigner::presolve_octets");
        //desres::profiler::Clock __(_);

        unsolved=_fragatoms;
        IdList keep;
        keep.reserve(_fragatoms.size());

        int npresolve=0;
        int nsolved=0;

        bool stillsolving;
        do {
            stillsolving=false;
            BOOST_FOREACH(Id aid1, unsolved){
                atom_t const& atm1=_mol->atom(aid1);
                int octval=DataForElement(atm1.atomic_number).maxOct/2;
                if(octval>4){
                    keep.push_back(aid1);
                    continue;
                } 
                int unkcount=0;
                electronMap::iterator lastatom=_atom_lp.end();
                electronMap::iterator lastbond=_bond_order.end();

                electronMap::iterator aiter=_atom_lp.find(aid1);
                if(aiter->second.lb==aiter->second.ub){
                    octval-=aiter->second.lb;
                }else{
                    unkcount++;
                    lastatom=aiter;
                } 
                IdList bonds=filteredBondsForAtom(_mol,aid1);
                BOOST_FOREACH(Id bid, bonds){
                    electronMap::iterator biter=_bond_order.find(bid);
                    if (biter->second.lb==biter->second.ub){
                        octval-=biter->second.lb;
                    }else{
                        unkcount++;
                        lastbond=biter;
                    }
                }
                if(unkcount==1){
                    // Get aid/bid that we can now solve for
                    if( lastatom != _atom_lp.end()){
                        electronRange &fixrange= lastatom->second;

                        if(!(octval >=fixrange.lb && octval <= fixrange.ub) && bonds.size()>0) {
                            std::stringstream msg;
                            msg << "Valence count solution is invalid for atom "<< lastatom->first <<
                                " ( " << octval << " " << fixrange.ub <<" )";
                            throw std::runtime_error(msg.str());
                        }

                        /* Only solve if possible lp count==0 or cant form hextet */
                        if(octval==0 || !allow_hextet_for_atom(aid1) ){
#if DEBUGPRINT
                            printf("AtomSolved: %d %d %d\n",lastatom->first,octval,fixrange.ub);
#endif
                            fixrange.lb=octval;
                            fixrange.ub=octval;
                            
                            stillsolving=true;
                            nsolved+=1;
                        }else{
                            keep.push_back(aid1);
                        }
                    }else{
                        electronRange &fixrange= lastbond->second;
#if DEBUGPRINT
                        bond_t const& _tmpbond=_mol->bond(lastbond->first);
                        printf("BondSolved: %d %d %d %d\n",_tmpbond.i,_tmpbond.j,octval,fixrange.ub);
#endif
                        if(!(octval >=fixrange.lb && octval <= fixrange.ub)) throw std::runtime_error("Bond order solution is invalid");
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
#if DEBUGPRINT 
        printf("Performed %d presolve iterations and finalized %d variables (%lu still unsolved)\n",npresolve,nsolved,unsolved.size());
#endif
    }

    ComponentAssigner::~ComponentAssigner(){
        lpsolve::delete_lp(_component_lp);
        lpsolve::delete_lp(_component_lpcopy);  
    }

     ComponentAssignerPtr ComponentAssigner::create(BondOrderAssignerPtr b, IdList const& comp, Id cid){
        //static desres::profiler::Symbol _("ComponentAssigner::initializeComponent");
        //desres::profiler::Clock __(_);

        ComponentAssignerPtr ca(new ComponentAssigner);

        assert(b!=NULL && comp.size()>0 );
        ca->_parent=b;
        ca->_component_atoms_present=std::set<Id>(comp.begin(), comp.end());
        ca->_component_id=cid;
        ca->_component_valence_count=0;
        // ca->_component_has_expanded_octets=false;
        ca->_component_solution_valid=false;

        return ca;
    }

    void ComponentAssigner::reset(){
        lpsolve::delete_lp(_component_lp);
        _component_lp=lpsolve::copy_lp(_component_lpcopy);
        _component_solution_valid=false;
    }

    void ComponentAssigner::setComponentCharge(int qTotal){
        if((_component_valence_count+qTotal)%2 != 0){
            std::stringstream ss;
            ss << "Desired charge of " <<qTotal<<" results in unpaired electrons in component "<<_component_id<<". ";
            ss << " (nelec="<<_component_valence_count+qTotal<<")\n";
            ss << "   Only closed shell configurations are supported";
            throw std::runtime_error(ss.str());
        }

        if(_component_solution_valid) reset();
        
        if(qTotal<0){
            lpsolve::set_bounds(_component_lp, _component_charge_col  , -qTotal,-qTotal);
            lpsolve::set_bounds(_component_lp, _component_charge_col+1, 0,      0);
            lpsolve::set_constr_type(_component_lp,_component_charge_row  , ROWTYPE_EQ);
            lpsolve::set_constr_type(_component_lp,_component_charge_row+1, ROWTYPE_LE);
        }else{
            lpsolve::set_bounds(_component_lp, _component_charge_col  , 0,      0);
            lpsolve::set_bounds(_component_lp, _component_charge_col+1, qTotal, qTotal);
            lpsolve::set_constr_type(_component_lp,_component_charge_row  , ROWTYPE_LE);
            lpsolve::set_constr_type(_component_lp,_component_charge_row+1, ROWTYPE_EQ);
        }
    }

    void ComponentAssigner::unsetComponentCharge(){
        if(_component_solution_valid) reset();
        BondOrderAssignerPtr parent=_parent.lock();
        lpsolve::set_bounds(_component_lp, _component_charge_col,  0, parent->max_component_charge);
        lpsolve::set_bounds(_component_lp, _component_charge_col+1,0, parent->max_component_charge);

        lpsolve::set_constr_type(_component_lp,_component_charge_row  , ROWTYPE_LE);
        lpsolve::set_constr_type(_component_lp,_component_charge_row+1, ROWTYPE_LE);
    }

    int ComponentAssigner::getSolvedComponentCharge(){
        if(!_component_solution_valid){
            std::stringstream msg;
            msg << "Cannot getSolvedComponentCharge from component "<<_component_id<<
                " with invalid integer linear program solution." <<
                " Did you call solveIntegerLinearProgram first?";
            throw std::runtime_error(msg.str());
        }
        int nrows=lpsolve::get_Norig_rows(_component_lp);
        int qTotal=
            - lpsolve::get_var_primalresult(_component_lp, nrows + _component_charge_col)
            + lpsolve::get_var_primalresult(_component_lp, nrows + _component_charge_col+1);
        return qTotal;
    }

    bool ComponentAssigner::solveComponentIntegerLinearProgram(){
        //static desres::profiler::Symbol _("ComponentAssigner::solveComponentIntegerLinearProgram");
        //desres::profiler::Clock __(_);

        if(_component_solution_valid) return _component_solution_valid;
        
        int status=lpsolve::solve(_component_lp);
        _component_solution_valid=(status==OPTIMAL || status==PRESOLVED);
#if DEBUGPRINT
        if(_component_solution_valid){
            int qTotal=getSolvedComponentCharge();
            printf("Solution was found for component %u with Charge= %d   objf= %6.3f\n",
                   _component_id, qTotal, lpsolve::get_objective(_component_lp));
            lpsolve::write_LP(_component_lp,stdout);
            lpsolve::print_objective(_component_lp);               
            lpsolve::print_solution(_component_lp,1);
        }else{
            printf("No solution found for component %u\n",_component_id);
            lpsolve::write_LP(_component_lp,stdout);
        }
#endif
        if(!_component_solution_valid) reset(); 
        
        return _component_solution_valid;
    }


    double ComponentAssigner::getSolvedComponentObjective(){
        if(!_component_solution_valid){
            std::stringstream msg;
            msg << "Cannot getSolvedComponentObjective from component "<<_component_id<< 
                " with invalid integer linear program solution." <<
                " Did you call solveIntegerLinearProgram first?";
            throw std::runtime_error(msg.str());
        }
        return lpsolve::get_objective(_component_lp);
    }

  
    void ComponentAssigner::check_resonance_path(resPathInfo const& respath, 
                                                 int direction, bool xferLonePair,
                                                 std::vector<int> const& soln,         
                                                 std::vector<int> & newsoln ){
        //static desres::profiler::Symbol _("ComponentAssigner::check_resonance_path");
        //desres::profiler::Clock __(_);
        static const int minbo=1;

        assert(abs(direction)==1);
        newsoln.clear();

        int idelta=direction*(xferLonePair ? 1 : -1);

        bool good=true;
        int deltab=idelta;
        /* quick check for valid resonance path */
        for (size_t idx=0; idx<respath.bondColsMaxOrders.size(); ++idx){
            std::pair<int,int> const& bdata=respath.bondColsMaxOrders[idx];
            int bcol  = bdata.first;
            int maxbo = bdata.second;

            int newbo=soln[bcol]+deltab;
            if(newbo<minbo || newbo>maxbo){
                good=false;
                break;
            }
            deltab*=-1;
        }
        
        if(!good) return;

        /* now do all the work */
        Id dcol,dqcol,acol,aqcol;
        if(direction==1){
            dcol    =respath.atom0_lpcol;
            dqcol   =respath.atom0_qcol;            
            acol    =respath.atom1_lpcol;
            aqcol   =respath.atom1_qcol; 
        }else{
            dcol    =respath.atom1_lpcol;
            dqcol   =respath.atom1_qcol;   
            acol    =respath.atom0_lpcol;
            aqcol   =respath.atom0_qcol; 
        }

#if DEBUGPRINT
        printf("found donor/acceptor path from %s to %s via:",
               lpsolve::get_origcol_name(_component_lp, dcol),
               lpsolve::get_origcol_name(_component_lp, acol));
#endif

        /* resonance path checks out, convert it to a new solution vector */
        newsoln=soln;
        if( xferLonePair ){
            // update valence electron pair counts
            newsoln[dcol]-=1;
            newsoln[acol]+=1;
        }

        // update bond orders
        deltab=idelta;
        for (size_t idx=0; idx<respath.bondColsMaxOrders.size(); ++idx){
            std::pair<int,int> const& bdata=respath.bondColsMaxOrders[idx];
            int bcol  = bdata.first;
            newsoln[bcol] += deltab;
            deltab *= -1;
        }
#if DEBUGPRINT
        for (size_t idx=0; idx<respath.bondColsMaxOrders.size(); ++idx){
            std::pair<int,int> const& bdata=respath.bondColsMaxOrders[idx];
            int bcol = bdata.first;
            printf("   %s", lpsolve::get_origcol_name(_component_lp, bcol));
        }
        printf("\n");
#endif
        
        /* update charge constraints
           NOTE: to avoid re-calculating the charge of the atom the assumed transitions are:
              donor  -1-> 0,  0->+1
           acceptor   0->-1, +1-> 0
           the charge constraints enters as n*q- and n*q+ (eg -1 charge is 1*q-). Under the above assumption, 
           we only need to perform binary flips of the charge constraints */
        if(newsoln[dqcol] > 1 || newsoln[dqcol+1]!=0){
            throw std::runtime_error("Bad charge constraint assumption. There is probably a bug in the code.");
        }else if(newsoln[dqcol]==1){
            newsoln[dqcol]=0;
        }else{
            newsoln[dqcol+1]=1;
        }
        if(newsoln[aqcol+1] > 1 || newsoln[aqcol]!=0){
            throw std::runtime_error("Bad charge constraint assumption. There is probably a bug in the code.");
        }else if(newsoln[aqcol+1]==1){
            newsoln[aqcol+1]=0;
        }else{
            newsoln[aqcol]=1;
        }

    }

    // return values are: antiaromatic = -1, non-aromatic = 0, aromatic = 1
    // For aromatic systems, also potentially returns newsoln with alternative bond orders.
    // NOTE: This routine only generates resonance forms for aromatic systems with alternating single/double bonds
    unsigned ComponentAssigner::check_aromatic_path(std::vector<ringAtomInfo> const& ringdata, 
                                                    std::vector<int> const& soln,
                                                    std::vector<int> & newsoln){

        //static desres::profiler::Symbol _("ComponentAssigner::check_aromatic_path");
        //desres::profiler::Clock __(_);

        newsoln.clear();

        /* Balaban, Chem Rev, 2004, 104, 2777-2812 
           Implementation of *simple* aromaticity detection
           Given ring atom types {X[2e], Y[1e], Z[0e]}, a ring of size M should satisfy:
           X+Y+Z=M
           X+0.5Y=2n+1  -> aromatic
           X+0.5Y=2n    -> antiaromatic
           The following additional restrictions are *NOT* considered here:
           1) sum of aromaticity constants should be -200< k[A] <200
           2) adjacent atoms of same type (for X and Z) destabilize the ring 
        */
        std::vector<unsigned> typecounts(AromaticAtom::INVALID);

        BOOST_FOREACH(ringAtomInfo const& atomInfo,ringdata){
            // Total Number of Bonds
            unsigned nb=atomInfo.nBonds;  
            // Number of lone pairs
            unsigned a0=soln[atomInfo.lonePairIdx]; 
            // Bond Order to previous ring atom
            unsigned b0=soln[atomInfo.bondIdxPrevious]; 
            // Bond Order to next ring atom
            unsigned b1=soln[atomInfo.bondIdxNext]; 
            // Bond Order from carbon to less electronegative exocyclic atom
            unsigned be=atomInfo.bondIdxExoCyclic==MIN_INVALID_ILP_COL ? 0 : soln[atomInfo.bondIdxExoCyclic];

            AromaticAtom::Type atype=AromaticAtom::Classify(nb,a0,b0,b1,be);
            if(atype>=AromaticAtom::INVALID){
                return AromaticRing::NONAROMATIC;
            }
            typecounts[atype]++;
        }

        AromaticRing::Type aro = AromaticRing::Classify(
                typecounts[AromaticAtom::X_TYPE],
                typecounts[AromaticAtom::Y_TYPE],
                typecounts[AromaticAtom::YEXT_TYPE],
                typecounts[AromaticAtom::Z_TYPE]);

        /* Try to generate new resonance structure if:
           1) ring is "aromatic", 
           2) has no internal lonepairs (nX==0),
           3) has no C=X external bonds(nYe==0), 
           4) no noncontributing atoms (nZ==0)
           This is equivalent to the below */
        if(aro==AromaticRing::AROMATIC && 
           typecounts[AromaticAtom::X_TYPE]==0 &&
           typecounts[AromaticAtom::YEXT_TYPE]==0 &&         
           typecounts[AromaticAtom::Z_TYPE]==0 ){
            int neworder=1; // flip flops between 1 & 2
            newsoln=soln;
            if(newsoln[ringdata[0].bondIdxNext]==1) neworder=2;
            BOOST_FOREACH(ringAtomInfo const& atomInfo, ringdata){
                /* Flip the bond order */
                newsoln[atomInfo.bondIdxNext]=3-newsoln[atomInfo.bondIdxNext];
                /* Check that new order is as expected */
                if(newsoln[atomInfo.bondIdxNext] != neworder){
                    newsoln.clear();
                    break;
                }
                /* Toggle expected order for next iteration */
                neworder=3-neworder;
            }
        }
        return aro;
    }


    
    void ComponentAssigner::initialize_donor_acceptors(std::vector<unsigned> & da_state,  
                                                      std::vector<int> & transfertype,
                                                      std::vector<resPathInfo> & resonance_paths ){
        //static desres::profiler::Symbol _("ComponentAssigner::initialize_donor_acceptors");
        //desres::profiler::Clock __(_);
        
        static const int _npairs=8;
        static const int _ndata=7;
        static const int da_data[_npairs][_ndata]= 
            {
                /* LD/LA-> lonepair donor/acceptor 
                   BD/BA-> bond donor/acceptor
                   type0 -> LD -> LA only, 
                   donorstates: LD=0, LA=1 
                   type1 -> LD -> LA,BD -> BA
                   donorstates: LD=2, LABD=3, BA=4
                   Data layout is: anum, #lp, #bonds, type, bo1...bo#bonds
                */                                              
                /*                                These have octet,              so do these,         but these have sextet */
                {6,  1, 3, 1, 1, 1, 1},   // LD: [  :C-](-R)(-R)(-R) -> LA,BD: [  C ](-R)(-R)(=R)  -> BA: [  C+](-R)(-R)(-R)
                {7,  1, 3, 0, 1, 1, 1},   // LD: [  :N ](-R)(-R)(-R) -> LA   : [  N+](-R)(-R)(=R)
                {7,  2, 2, 1, 1, 1, 0},   // LD: [ ::N-](-R)(-R)     -> LA,BD: [ :N ](-R)(=R)      -> BA: [ :N+](-R)(-R)
                {7,  2, 1, 1, 2, 0, 0},   // LD: [ ::N-](=R)         -> LA,BD: [ :N ](#R)          -> BA: [ :N+](=R)
                {8,  2, 2, 0, 1, 1, 0},   // LD: [ ::O ](-R)(-R)     -> LA   : [  O+](-R)(=R)
                {8,  3, 1, 0, 1, 0, 0},   // LD: [:::O-](-R)         -> LA,BD: [::O ](=R)          {SKIP-> BA: [::O+](-R)}
                {16, 2, 2, 0, 1, 1, 0},   // LD: [ ::S ](-R)(-R)     -> LA   : [ :S+](-R)(=R) 
                {16, 3, 1, 0, 1, 0, 0}    // LD: [:::S-](-R)         -> LA,BD: [::S ](=R)          {SKIP-> BA: [::S+](-R)}
            };
        
        da_state.clear();
        transfertype.clear();
        resonance_paths.clear();

        bool allow_hextet_resonance=false;

        BondOrderAssignerPtr parent=_parent.lock();

        IdList da_atoms;
        BOOST_FOREACH(ilpMap::value_type const& apair, _component_atom_cols){
            Id aid=apair.first;
            int atomcol=apair.second;
            std::vector<int> orders;
            atom_t const& ai=parent->_mol->atom(aid);
            for (int idx=0; idx<_npairs; ++idx){
                int anum=da_data[idx][0];
                int nlp=da_data[idx][1];
                Id nbonds=da_data[idx][2];
                IdList bonds=filteredBondsForAtom(parent->_mol,aid);
                if (!(ai.atomic_number == anum && bonds.size() == nbonds )) continue;
                int datype=da_data[idx][3];
                orders.clear();

                BOOST_FOREACH(Id bid, bonds){
                    orders.push_back(_component_solution[asserted_find(_component_bond_cols,bid)]);
                }
                std::sort(orders.begin(),orders.end());
                bool boequal=std::equal(orders.begin(),orders.end(),&da_data[idx][4]);
                bool istype1=datype && allow_hextet_resonance;
                if(_component_solution[atomcol] == nlp && boequal ){
                    da_atoms.push_back(aid);
                    da_state.push_back(istype1 ? 2 : 0);
                }else if(_component_solution[atomcol] == nlp-1){
                    if(istype1 && boequal){
                        da_atoms.push_back(aid);
                        da_state.push_back(4);
                    }else{
                        orders[orders.size()-1]--;
                        if(std::equal(orders.begin(),orders.end(),&da_data[idx][4]) ){
                            da_atoms.push_back(aid);
                            da_state.push_back(istype1 ? 3 : 1);
                        }
                    }
                }
            }
        }
        /* generate all possible paths from potential "donor"/"acceptor" atom to "acceptor"/"donor" atom. 
           Only generated in one direction as the reverse path can be used when necessary 
         */
        unsigned dacount=da_atoms.size();
        for (unsigned i=0; i<dacount; ++i){
            Id donor=da_atoms[i];
            for (unsigned j=i+1; j<dacount; ++j){
                Id acceptor=da_atoms[j];

                /* An electron transfer between disimilar atoms will change 
                   the objective so we dont need to check them out 
                   (we are only interested in equivalent "best" solutions */
                if(parent->_mol->atom(donor).atomic_number != 
                   parent->_mol->atom(acceptor).atomic_number) continue;

                std::set<Id> visited;
                typedef std::map<Id, Id> pathMap;
                pathMap path;
                
                typedef std::pair<Id, Id> node;
                std::deque<node> stack;
                
                stack.push_back(node(donor,BadId));
                while(! stack.empty() ){
                    node &current_node=stack.back();
                    Id curatom=current_node.first;
                    Id prevatom=current_node.second;            
                    stack.pop_back();
                    
                    pathMap::iterator piter=path.find(prevatom);
                    // Unroll visited and current path to correct state
                    while( piter != path.end() ){
                        Id next=piter->second;
                        visited.erase(next);
                        path.erase(piter);
                        piter=path.find(next);
                    }
                    
                    path.insert(pathMap::value_type(prevatom,curatom));
                    visited.insert(curatom);
                    
                    if(curatom==acceptor){
                        IdList apath;
                        // traverse path (sequence of visited atoms)
                        piter=path.find(donor);
                        while(piter!=path.end()){
                            apath.push_back(piter->first);
                            piter=path.find(piter->second);
                        }
                        apath.push_back(acceptor);
                        assert(apath.front()==donor && apath.back()==acceptor);
                        // path length (# of atoms in path) must be odd for valid electron transfer
                        if(apath.size()%2 == 1){
#if DEBUGPRINT2
                            printf("Possible Atom Resonance Path from %u to %u\n",donor,acceptor);
#endif
                            resonance_paths.push_back(resPathInfo());
                            resPathInfo &resPath=resonance_paths.back();

                            resPath.atom0_daidx=i;
                            resPath.atom0_lpcol=asserted_find(_component_atom_cols, donor);   
                            resPath.atom0_qcol =asserted_find(_component_atom_charge_cols, donor);
                   
                            resPath.atom1_daidx=j;
                            resPath.atom1_lpcol=asserted_find(_component_atom_cols, acceptor);   
                            resPath.atom1_qcol =asserted_find(_component_atom_charge_cols, acceptor);  
                               
                            for( size_t idx=0;idx<apath.size()-1;++idx){
                                resPath.bondColsMaxOrders.push_back(std::pair<int,int>());
                                std::pair<int,int> &bondData=resPath.bondColsMaxOrders.back();
                                Id bid=parent->_mol->findBond(apath[idx],apath[idx+1]);
                                bondData.first=asserted_find(_component_bond_cols,bid);
                                int maxbo=parent->max_bond_order(apath[idx]);
                                int maxbo2=parent->max_bond_order(apath[idx+1]);
                                if(maxbo>maxbo2) maxbo=maxbo2;
                                bondData.second=maxbo;
                            }
                        }
                    }else{
                        IdList bonded=filteredBondedAtoms(parent->_mol, curatom);
                        BOOST_FOREACH( Id nextatom, bonded){
                            /* A resonance path can only exist where the initial bond orders were indeterminate */
                            if(visited.count(nextatom) > 0 || 
                               _component_atom_cols.find(nextatom) == _component_atom_cols.end()) continue;
                            
                            stack.push_back(node());
                            node &next_node=stack.back();
                            next_node.first=nextatom;
                            next_node.second=curatom;
                        }
                    }
                }
            }
        }
#if DEBUGPRINT
        printf("Found %zu resonance path candidates\n", resonance_paths.size());
#endif
        
        transfertype.reserve(25);
        for (int i=0;i<5;++i){
            for(int j=0;j<5;++j){
                if ( ((i==0 || i==2) && (j==1 || j==3)) || (i==3 && j==4) ) {
                    // Atom i donor, atom j acceptor
                    transfertype.push_back(1);
                }else if ( ((j==0 || j==2) && (i==1 || i==3)) || (j==3 && i==4) ) {
                    // Atom j donor, atom i acceptor
                    transfertype.push_back(-1);
                }else{
                    // Invalid donor/acceptor combination
                    transfertype.push_back(0);          
                }
            }
        }
    }

    void ComponentAssigner::initialize_aromatic_ringdata(std::vector< std::vector< ringAtomInfo > > &rings_to_process){
        //static desres::profiler::Symbol _("ComponentAssigner::initialize_aromatic_ringdata");
        //desres::profiler::Clock __(_);

        BondOrderAssignerPtr parent=_parent.lock();
        rings_to_process.clear();
        /* Initialize rings for resonance form generator */
        BOOST_FOREACH(ilpMap::value_type const& rpair, _component_ring_cols){
            Id rid=rpair.first;
            size_t natoms= parent->_planar_rings[rid].second.size()-1;
            Id previous= parent->_planar_rings[rid].second[natoms-1];

            rings_to_process.push_back(std::vector<ringAtomInfo>());
            std::vector<ringAtomInfo> &ringdata=rings_to_process.back();
            ringdata.reserve(natoms);
            for(size_t iatom=0; iatom<natoms; ++iatom){
                ringdata.push_back(ringAtomInfo());
                ringAtomInfo &atomInfo=ringdata.back();

                Id current = parent->_planar_rings[rid].second[iatom];
                Id next =  parent->_planar_rings[rid].second[iatom+1];

                IdList bonds=filteredBondsForAtom(parent->_mol,current);
                Id nbonds=bonds.size();
                assert(nbonds<4);
                atomInfo.aid=current;
                atomInfo.nBonds=nbonds;
                atomInfo.lonePairIdx=asserted_find( _component_atom_cols, current);

                int prevcol=MIN_INVALID_ILP_COL;
                int nextcol=MIN_INVALID_ILP_COL;
                int exocol=MIN_INVALID_ILP_COL;
                BOOST_FOREACH( Id bid, bonds){
                    Id other=parent->_mol->bond(bid).other(current);
                    if(other==previous){
                        prevcol=asserted_find( _component_bond_cols, bid);
                    }else if(other==next){
                        nextcol=asserted_find( _component_bond_cols, bid);
                    }else if(parent->_mol->atom(current).atomic_number==6 && 
                             parent->_mol->atom(other).atomic_number == 6){
                        exocol=asserted_find( _component_bond_cols, bid);
                    }
                }
                assert(nextcol!=MIN_INVALID_ILP_COL && prevcol!=MIN_INVALID_ILP_COL);
                atomInfo.bondIdxPrevious=prevcol;
                atomInfo.bondIdxNext=nextcol;
                atomInfo.bondIdxExoCyclic=exocol;    

                previous=current;
            }
        }
    }

    void ComponentAssigner::update_aring_constraint_for_resonance_form(std::vector<std::vector<double> > const& arom_cons,
                                                                       std::vector<int> const& old_soln, 
                                                                       std::vector<int> & new_soln){
        //static desres::profiler::Symbol _("ComponentAssigner::update_aring_constraint_for_resonance_form");
        //desres::profiler::Clock __(_);

        assert(new_soln.size()==old_soln.size());
        assert(arom_cons.size()==_component_ring_cols.size());

        size_t ncols=old_soln.size();
        size_t idx=0;
        BOOST_FOREACH(ilpMap::value_type const& rpair, _component_ring_cols){
            std::vector<double> const& vals=arom_cons[idx];
            assert(ncols==vals.size());
            double oldsum=0;
            double newsum=0;
            for (size_t vidx=1;vidx<vals.size();++vidx){
                if(vals[vidx]==0.0) continue;
                oldsum+=old_soln[vidx]*vals[vidx];
                newsum+=new_soln[vidx]*vals[vidx];
            }
            /* number of electrons gained/lost in this ring versus initial solution */
            int delta=double_to_int(newsum-oldsum,0.0);

            int colid=rpair.second;
            int nearest_n=old_soln[colid];
            /*  desum is the electron change needed to be aromatic */
            int desum=old_soln[colid+2]-old_soln[colid+1];
            /* *subtract* delta from desum to get new additional ecount needed to be aromatic
               eg, positive 'desum'=needs electrons to be aromatic, positive 'delta'=gained electron 
            */
            desum-=delta;
            /* bring desum within range by adjusting nearest_n */
            while(desum<-2){
                nearest_n+=1;
                desum+=4;
            }
            while(desum>2){
                nearest_n-=1;
                desum-=4;
            }
            /* nearest_n cannot be negative */
            assert(nearest_n>=0);

            /* update solution vector with corrected ring variables */
            if(desum<0){
                new_soln[colid  ]=nearest_n;
                new_soln[colid+1]=-desum;
                new_soln[colid+2]=0;                
            }else if (desum>0){
                new_soln[colid  ]=nearest_n;
                new_soln[colid+1]=0;
                new_soln[colid+2]=desum;  
            }else{
                new_soln[colid  ]=nearest_n;
                new_soln[colid+1]=0;
                new_soln[colid+2]=0;
            }
            ++idx;
        }
    }

    void ComponentAssigner::generate_resonance_forms_to_check(std::vector<std::vector<double> > const& arom_cons,
                                                              std::set< std::vector<int> > &rforms){
        //static desres::profiler::Symbol _("ComponentAssigner::generate_resonance_forms_to_check");
        //desres::profiler::Clock __(_);

        std::vector<unsigned> da_state;
        std::vector<int> transferDirection;
        std::vector<resPathInfo> resonancePaths;
        initialize_donor_acceptors(da_state, transferDirection, resonancePaths );

        std::vector< std::vector< ringAtomInfo > >ringPaths;
        initialize_aromatic_ringdata(ringPaths);

        typedef std::map< std::vector<int>, std::vector<unsigned> > resformMap;
        resformMap rforms_to_process;
        resformMap::iterator rform,rformtmp;
        rformtmp=rforms_to_process.insert(resformMap::value_type(_component_solution,resformMap::mapped_type())).first;
        rformtmp->second.swap(da_state);

        int nforms=1;
        std::vector<int> params;
        std::vector<int> newparams;

        rforms.clear();
        while(!rforms_to_process.empty()){
#if DEBUGPRINT
            printf("Processing new primary rform\n");
#endif
            // pull first value from map
            rform=rforms_to_process.begin();
            params=rform->first;
            da_state.swap(rform->second);
            rforms_to_process.erase(rform);
            // Add to seen resonance forms
            rforms.insert(params);

            // Find electron resonance paths from donor to acceptors
            newparams.clear();
            BOOST_FOREACH(resPathInfo const& respath, resonancePaths){
                unsigned state0=da_state[respath.atom0_daidx];
                unsigned state1=da_state[respath.atom1_daidx];
                int direction=transferDirection[5*state0+state1];
                if(direction==0) continue;
                bool xferlp=!(state0>2 && state1>2);

                check_resonance_path(respath, direction, xferlp, params, newparams);
                if(newparams.size()>0){
                    update_aring_constraint_for_resonance_form(arom_cons, params, newparams);
                    if(rforms.find(newparams)==rforms.end()){
                        rformtmp=rforms_to_process.lower_bound(newparams);
                        if(rformtmp==rforms_to_process.end() || rforms_to_process.key_comp()(newparams,rformtmp->first)){
                            nforms++;
                            
                            /* fill vector after insertion */
                            rformtmp=rforms_to_process.insert(rformtmp, resformMap::value_type(newparams,resformMap::mapped_type()));
                            rformtmp->second=da_state;
                            if(direction==1){
                                rformtmp->second[respath.atom0_daidx]= state0+1;
                                rformtmp->second[respath.atom1_daidx]= state1-1;
                            }else{
                                rformtmp->second[respath.atom0_daidx]= state0-1;
                                rformtmp->second[respath.atom1_daidx]= state1+1;
                            }
                            
                        }  
                    }            
                }
            }

            // Find aromatic ring resonance paths (does not modify donor/acceptor)
            BOOST_FOREACH(std::vector<ringAtomInfo> const& ringdata, ringPaths){
                std::vector<int> newparams;

                check_aromatic_path(ringdata, params, newparams);
                if(newparams.size()>0){
                    update_aring_constraint_for_resonance_form(arom_cons, params, newparams);
                    if(rforms.find(newparams)==rforms.end() ){
                        rformtmp=rforms_to_process.lower_bound(newparams);
                        if(rformtmp==rforms_to_process.end() || rforms_to_process.key_comp()(newparams,rformtmp->first)){
                            nforms++;
                            
                            rformtmp=rforms_to_process.insert(rformtmp,resformMap::value_type(newparams,resformMap::mapped_type()));
                            rformtmp->second=da_state;
                            
                        }           
                    }
                }
            }
        }
#if DEBUGPRINT
        printf("# Resonance Forms Generated: %d\n",nforms);
        printf("# Resonance Forms Kept Initial: %ld\n",rforms.size());
#endif

    }
    
    void ComponentAssigner::sum_resonance_solutions(std::set< std::vector<int> > const& rforms){
        unsigned nvars=_component_objf.size();

        Quadsum best_obj=0.0;
        for (unsigned i=1; i<nvars;++i){
            best_obj += _component_solution[i]*_component_objf[i];
        } 

        std::vector< int > cur_soln(nvars,0);
        std::set< std::vector<int> >::const_iterator rforms_iter;
        for(rforms_iter=rforms.begin(); rforms_iter!= rforms.end(); ++rforms_iter){
            assert(rforms_iter->size()==nvars);

            Quadsum cur_obj=0;
            for (unsigned i=1; i<nvars;++i){
                cur_obj += (*rforms_iter)[i]*_component_objf[i];
            }

            double objdiff=cur_obj.result()-best_obj.result();

            if(objdiff>0.0){
                continue;
            }else if (objdiff<0.0){
                printf("Warning: resonance generator detected better solution... Something is probably wrong in solution vector\n");
                continue;
            }
            cur_soln[0]++;
            for (unsigned i=1; i<nvars;++i){
                cur_soln[i]+=(*rforms_iter)[i];
            }
        }

        double nforms=cur_soln[0] > 0.0 ?  cur_soln[0] : 1.0;
        _component_resonant_solution.resize(nvars);
        _component_resonant_solution[0]=nforms;  
        if(nforms>0.0){
            for (unsigned i=1; i<nvars;++i){
                _component_resonant_solution[i]=cur_soln[i]/nforms;
            }
        }else{
            for (unsigned i=1; i<nvars;++i){
                _component_resonant_solution[i]=_component_solution[i];
            }
        }

#if DEBUGPRINT
        printf("# Resonance Forms Kept Final: %d\n",static_cast<int>(nforms));
        printf("Final Params:        ");
        for (unsigned i=1; i<nvars;++i) printf("%3.2f*%.2f ",_component_resonant_solution[i],_component_objf[i]);
        printf("\n");
#endif
    }


    /* Helper function for extractComponentSolution */
    void ComponentAssigner::extract(ilpMap const& im, bool isCharge, solutionMap & sm){
        BOOST_FOREACH(ilpMap::value_type const& ipair, im ){       
            solutionMap::iterator entry=sm.lower_bound(ipair.first); 
            int nonres=_component_solution.at(ipair.second); 
            double res=_component_resonant_solution.at(ipair.second);
            if(isCharge){
                /* Total charge takes up 2 columns... one for + charge and one for - */
                nonres=_component_solution.at(ipair.second+1)-nonres;
                res=_component_resonant_solution.at(ipair.second+1)-res;
            }
            if(entry==sm.end() || sm.key_comp()(ipair.first,entry->first)){ 
                sm.insert(entry, solutionMap::value_type(ipair.first,solutionValues(nonres,res))); 
            }else{                                                 
                assert(nonres==entry->second.nonresonant && res==entry->second.resonant); 
            }                                                       
        }                        
    }


    void ComponentAssigner::get_ilp_solution(lpsolve::_lprec *lp, std::vector<int> &solution){

        int nrows=lpsolve::get_Norig_rows(lp);
        int ncols=lpsolve::get_Norig_columns(lp);
        unsigned solsize=ncols+1;

        solution.assign(solsize,0);
        for (unsigned idx=1; idx<solsize;++idx){
            double result=lpsolve::get_var_primalresult(lp, nrows + idx);
            solution[idx]=double_to_int(result, 0.0);
        }
    }


    void ComponentAssigner::extractComponentSolution(solutionMap &atominfo,
                                                     solutionMap &bondinfo,
                                                     solutionMap &chargeinfo){
        //static desres::profiler::Symbol _("ComponentAssigner::extractComponentSolution");
        //desres::profiler::Clock __(_);

        assert(_component_solution_valid);

        ComponentAssigner::get_ilp_solution(_component_lp,_component_solution);

        std::vector<std::vector<double> > arom_cons;
        BOOST_FOREACH(ilpMap::value_type const& rpair, _component_ring_cols){
            Id ridx=rpair.first;
            int colid =rpair.second;
            arom_cons.push_back(std::vector<double>());
            std::vector<double> &rowdata=arom_cons.back();
            double target;
            generate_ring_constraint(ridx, colid, target, rowdata);
        }
        
        std::set< std::vector<int> > rforms;
        generate_resonance_forms_to_check( arom_cons, rforms);
        sum_resonance_solutions( rforms );
        
        extract(_component_atom_cols, false, atominfo);
        extract(_component_bond_cols, false, bondinfo);
        extract(_component_atom_charge_cols, true, chargeinfo);
    }


    /* Helper function for adding a column to the integer linear program
       Sets bounds and objective function contribution (as well as column 
       label for debugging) */
    int ComponentAssigner::add_column_to_ilp(lpsolve::_lprec *lp, std::string const& colname, double penalty, double lb, double ub){
        double coldata[1]={penalty};
        lpsolve::add_column(lp,coldata);
        int colid=lpsolve::get_Norig_columns(lp);
        assert(colid>MIN_INVALID_ILP_COL);
        lpsolve::set_bounds(lp,colid,lb,ub);
        lpsolve::set_int(lp,colid,true);
        lpsolve::set_col_name(lp,colid,const_cast<char*>(colname.c_str()));
        return colid;
    }


    /* prefer lone pairs on more electronegative atom */
    void ComponentAssigner::set_atom_lonepair_penalties(){
        //static desres::profiler::Symbol _("ComponentAssigner::set_atom_penalties");
        //desres::profiler::Clock __(_);

        _component_atom_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();
        std::ostringstream ss;
        BOOST_FOREACH(Id aid1, _component_atoms_present){
            atom_t const& atm1=parent->_mol->atom(aid1);
            int anum1=atm1.atomic_number;
            electronRange const& range= asserted_find(parent->_atom_lp,aid1);
            std::ostringstream ss;
            ss.str("");
            ss << "a_"<<aid1;
            /* factor of 2 since there are 2 electrons per lp */
            double objv=-1*DataForElement(anum1).eneg;
            int colid=add_column_to_ilp(_component_lp,ss.str(),objv, range.lb, range.ub);
            _component_atom_cols.insert(std::make_pair(aid1,colid)); 
            
        }
    }

    /* prefer bonds between more electronegative atom pairs */
    void ComponentAssigner::set_bond_penalties(){
        //static desres::profiler::Symbol _("ComponentAssigner::set_bond_penalties");
        //desres::profiler::Clock __(_);

        _component_bond_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();
        /* partition factor must be in range */
        assert(parent->bond_eneg_factor>=0 && parent->bond_eneg_factor <=1);
        double large_factor=1.0-parent->bond_eneg_factor;

        std::ostringstream ss;
        BOOST_FOREACH(Id aid1, _component_atoms_present){
            int anum1=parent->_mol->atom(aid1).atomic_number;
            IdList bonds=filteredBondsForAtom(parent->_mol,aid1);
            BOOST_FOREACH(Id bid, bonds){
                electronRange const& range= asserted_find(parent->_bond_order,bid);
                // Have we already added this bond? 
                if ( _component_bond_cols.find(bid) != _component_bond_cols.end()) continue;
                Id aid2=parent->_mol->bond(bid).other(aid1);
                atom_t const& atm2=parent->_mol->atom(aid2);
                int anum2=atm2.atomic_number;
                ss.str("");
                ss << "b_"<<aid1<<"_"<<aid2;
                double objv;
                /* Total factor of 2 since there are 2 electrons per bond */
                if(DataForElement(anum1).eneg <= DataForElement(anum2).eneg){
                    objv=-2.0*(large_factor*DataForElement(anum2).eneg +
                               parent->bond_eneg_factor*DataForElement(anum1).eneg);
                }else{
                    objv=-2.0*(large_factor*DataForElement(anum1).eneg + 
                               parent->bond_eneg_factor*DataForElement(anum2).eneg);
                }
                int colid=add_column_to_ilp(_component_lp,ss.str(),objv,range.lb,range.ub);
                _component_bond_cols.insert(std::make_pair(bid,colid));
                
            }
        }
    }


    bool BondOrderAssigner::allow_hextet_for_atom(Id aid1){
        int anum1=_mol->atom(aid1).atomic_number;
        Id nbonds=filteredBondCountForAtom(_mol, aid1);
        /* Allow hextets for unsaturated carbon, nitrogen and oxygen */
        if ( (anum1==6 && (nbonds==3 || nbonds==2)) || 
             (anum1==7 && (nbonds==2 || nbonds==1)) ||
             (anum1==8 && (nbonds==1             ))) return true;
        return false;
    }

    /* penalize atoms for having a hextet */
    void ComponentAssigner::set_atom_hextet_penalties(){
        //static desres::profiler::Symbol _("ComponentAssigner::set_atom_hextet_penalties");
        //desres::profiler::Clock __(_);

        _component_atom_hextet_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();

        std::ostringstream ss;
        BOOST_FOREACH(Id aid1, _component_atoms_present){

            if ( !parent->allow_hextet_for_atom(aid1) ) continue;
            ss.str("");
            ss << "hex_"<<aid1;
            int colid=add_column_to_ilp(_component_lp,ss.str(),parent->atom_hextet_penalty,0,1);

            _component_atom_hextet_cols.insert(std::make_pair(aid1,colid));
      
        }
    }


    bool ComponentAssigner::gen_charge_penalty_for_atom(Id aid1){
        BondOrderAssignerPtr parent=_parent.lock();
        int anum1=parent->_mol->atom(aid1).atomic_number;
        Id nbonds=filteredBondCountForAtom(parent->_mol, aid1);
        /* no charge penalty necessary for ions, hydrogen, saturated carbons/nitrogen
           since charge is predetermined */
        if ( nbonds==0 || anum1==1 || (nbonds==4 && (anum1==6 || anum1==7)) ) return false;
        return true;
    }

    /* penalize atoms for having a charge */
    void ComponentAssigner::set_atom_charge_penalties(){
        //static desres::profiler::Symbol _("ComponentAssigner::set_atom_charge_penalties");
        //desres::profiler::Clock __(_);

        _component_atom_charge_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();

        /* offset is used to remap electronegativity scale into
           penalty for creating negativly charged atoms.
           TODO: Investigate this value */
        static const double offset=DataForElement(10).eneg+DataForElement(87).eneg;
 
        std::ostringstream ss;
        BOOST_FOREACH(Id aid1, _component_atoms_present){

            if ( !gen_charge_penalty_for_atom(aid1) ) continue;
        
            int anum1=parent->_mol->atom(aid1).atomic_number;
            ss.str("");
            ss << "qm_"<<aid1;
            /* smaller penalty to put '-' charge on more electronegative atom */
            double objv= offset-DataForElement(anum1).eneg;
            int colid=add_column_to_ilp(_component_lp,ss.str(),parent->atom_charge_penalty_factor*objv,0,parent->max_atom_charge);
            
            /* We only need to keep track of the first column id */
            _component_atom_charge_cols.insert(std::make_pair(aid1,colid));
            
            ss.str("");
            ss << "qp_"<<aid1;
            /* smaller penalty to put '+' charge on less electronegative atom */
            objv=DataForElement(anum1).eneg;
            add_column_to_ilp(_component_lp,ss.str(),parent->atom_charge_penalty_factor*objv,0,parent->max_atom_charge);
        }
    }

    /* penalize rings for not being aromatic */
    void ComponentAssigner::set_aromatic_ring_penalties(){
        
        _component_ring_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();

        /* penalty proportional to carbon eneg for not forming aromatic ring
           TODO: Investigate value */
        static const double penalty=DataForElement(6).eneg;
        double objv=parent->ring_penalty_factor*penalty;

        std::ostringstream ss;
        // Try to make planar rings aromatic 
        for(Id ridx=0; ridx<parent->_planar_rings.size(); ++ridx){
            IdList const& ring=parent->_planar_rings[ridx].second;
            bool addRing=true;
            BOOST_FOREACH(Id aid1, ring){
                if(filteredBondCountForAtom(parent->_mol,aid1) >3 || _component_atoms_present.count(aid1)==0){
                    addRing=false;
                    break;
                }
            }
            if(!addRing) continue;

            ss.str("");
            ss << "rn_"<<ridx; // 'n' in 4n+2 pi electrons for aromaticity
            int colid=add_column_to_ilp(_component_lp,ss.str(),0,0,1000);

            /* Only need to keep track of the first column id */
            _component_ring_cols.insert(std::make_pair(ridx,colid));

            ss.str("");
            ss << "rs_"<<ridx; // do we need to subtract electrons from ring to make it aromatic?
            add_column_to_ilp(_component_lp,ss.str(), objv, 0, 2);

            ss.str("");
            ss << "ra_"<<ridx; // or add electrons to ring to make it aromatic?
            add_column_to_ilp(_component_lp,ss.str(), objv, 0, 2);

        }
    }

    /* penalize components for having a charge */
    void ComponentAssigner::set_component_charge_penalty(){

        BondOrderAssignerPtr parent=_parent.lock();

        std::ostringstream ss;
        double maxeneg=0.0;
        double mineneg=DataForElement(10).eneg;
        BOOST_FOREACH(Id  aid1, _component_atoms_present){
            double eneg = DataForElement(parent->_mol->atom(aid1).atomic_number).eneg;
            if(eneg < mineneg) mineneg=eneg;
            if(eneg > maxeneg) maxeneg=eneg;
        }
        ss.str("");
        ss << "qCompm_"<< _component_id; 
        /* assume negative charge gets added to atom with largest electronegativity */
        int colid=add_column_to_ilp(_component_lp,ss.str(), parent->component_charge_penalty_factor*maxeneg,0,parent->max_component_charge);
        
        /* Only need to keep track of the first column id */
        _component_charge_col=colid;
        
        ss.str("");
        ss << "qCompp_"<< _component_id; 
        /* assume electron gets removed from atom with smallest electronegativity */
        add_column_to_ilp(_component_lp,ss.str(),parent->component_charge_penalty_factor*mineneg,0,parent->max_component_charge);
        
    }

    void ComponentAssigner::break_ring_symmetry(){

        BondOrderAssignerPtr parent=_parent.lock();

        /* Add a *small* purtubation to the ringsize/2 shortest bonds in rings 
           with identical bonds. If the purturbation is too large,
           it could adversly impact the initial solution */
        static const double scale=1.0001;

        BondOrderAssigner::RingList* rlists[2]={&parent->_planar_rings,&parent->_nonplanar_rings};

        for (size_t ridx=0; ridx<2; ++ridx){
            BOOST_FOREACH(BondOrderAssigner::RingPair const& rp, *rlists[ridx]){
                IdList const& ring=rp.second;
                // (first and last atom in rings are the same... (size-1)
                size_t ringsize=ring.size()-1;
                bool purturbRing=true;
                int targetanum=parent->_mol->atom(ring[0]).atomic_number;
                BOOST_FOREACH(Id  aid1, ring){
                    if(filteredBondCountForAtom(parent->_mol,aid1) >3 || 
                       _component_atoms_present.count(aid1)==0 ||
                       parent->_mol->atom(aid1).atomic_number != targetanum){
                        purturbRing=false;
                        break;
                    }
                }
                if(!purturbRing) continue;
                
                std::multimap<float,Id> distanceMap;
                for(IdList::const_iterator miter=ring.begin(); miter != ring.end()-1; ++miter){
                    Id current = *miter;
                    Id next = *(miter+1);
                    atom_t const& atm1=parent->_mol->atom(current);
                    atom_t const& atm2=parent->_mol->atom(next);
                    double dx=atm2.x - atm1.x;
                    double dy=atm2.y - atm1.y;
                    double dz=atm2.z - atm1.z;
                    float dist=sqrt(dx*dx+dy*dy+dz*dz);
                    Id bid =parent->_mol->findBond(current, next);
                    distanceMap.insert(std::make_pair(dist,bid));
                }
                std::multimap<float,Id>::const_iterator dmiter=distanceMap.begin();
                ringsize/=2;
                for(size_t idx=0; idx<=ringsize; ++idx,++dmiter){
                    int colid=asserted_find(_component_bond_cols,dmiter->second);
                    double valold=lpsolve::get_mat(_component_lp,0,colid);
                    lpsolve::set_mat(_component_lp,0,colid, scale*valold);
                }            
            }
        }
    }

    
    void ComponentAssigner::generate_ring_constraint(Id ridx, int rcolid, double &target,
                                                     std::vector<double> &rowdata){
        BondOrderAssignerPtr parent=_parent.lock();

        int ncols=lpsolve::get_Norig_columns(_component_lp);
        rowdata.assign(ncols+1,0);
        
        target= 2.0; // pi electrons, (4*n + 2) electrons
        rowdata.at(rcolid)=-4;
        
        IdList const& ring=parent->_planar_rings.at(ridx).second;
        size_t ringsize=ring.size()-1;
        Id previous= ring[ringsize-1];
        IdList::const_iterator miter;
        for(miter=ring.begin(); miter != ring.end()-1; ++miter){
            Id current = *miter;
            IdList bonds=filteredBondsForAtom(parent->_mol, current);
            Id nbonds=bonds.size();
            assert(nbonds<4);
            
            /* two electrons per lone pair... */
            rowdata.at(asserted_find(_component_atom_cols,current))=2;
            /* corrected for how many are available in the pi system */
            target+=2*(3 - nbonds);
            
            Id next = *(miter+1);
            BOOST_FOREACH( Id bid, bonds){
                Id other=parent->_mol->bond(bid).other(current);
                if(other==previous){
                    /* Nothing to do for other==previous */
                }else if(other==next){
                    /* two electrons per bond order in ring... */
                    rowdata.at(asserted_find(_component_bond_cols,bid))=2;
                    /* corrected for how many are available in the pi system */
                    target+=2;
                }else if(parent->_mol->atom(current).atomic_number==6 && 
                         DataForElement(6).eneg >= DataForElement(parent->_mol->atom(other).atomic_number).eneg){
                    /* endocyclic carbon to exocyclic atom (could be double bond)
                       only *one* unshared electron per bond order
                       need += because 'exocyclic' bond may actually be an internal bond (naphthlene)
                    */ 
                    rowdata.at(asserted_find(_component_bond_cols, bid))+=1;
                    /* corrected by how many are available to the pi system */
                    target+=1;        
                }
            }
            previous=current;
        }
    }

    void ComponentAssigner::add_aromatic_ring_constraints(){

        int ncols=lpsolve::get_Norig_columns(_component_lp);
        double target;
        std::vector<double> rowdata;

        /* colid   : n in 4n+2 pi electrons
           colid+1 : electrons need to be removed to make ring aromatic
           colid+2 : electrons need to be added   to make ring aromatic
        */
        BOOST_FOREACH(ilpMap::value_type const& rpair, _component_ring_cols){
            Id ridx=rpair.first;
            int colid =rpair.second;
            generate_ring_constraint(ridx, colid, target, rowdata);
            
            /* constraint on subtracting electrons from ring */
            rowdata.at(colid+1)=-1;
            rowdata.at(colid+2)= 0;            
            lpsolve::add_constraint(_component_lp, &rowdata[0], ROWTYPE_LE, target);

            /* constraint on adding electrons to ring */
            for(int k=0;k<ncols;++k) rowdata.at(k+1)*=-1;
            rowdata.at(colid+1)= 0;
            rowdata.at(colid+2)=-1;     
            lpsolve::add_constraint(_component_lp, &rowdata[0], ROWTYPE_LE,-target);
        }

    }


    void ComponentAssigner::add_atom_octet_and_charge_constraints(){
        //static desres::profiler::Symbol _("ComponentAssigner::add_atom_octet_and_charge_constraints");
        //desres::profiler::Clock __(_);

        BondOrderAssignerPtr parent=_parent.lock();

        int ncols=lpsolve::get_Norig_columns(_component_lp);
        std::vector<double> rowdata;

        BOOST_FOREACH(ilpMap::value_type const& apair, _component_atom_cols){
            rowdata.assign(ncols+1,0);
            Id aid1 = apair.first;
            IdList bonds=filteredBondsForAtom(parent->_mol, aid1);
            Id nbonds=bonds.size();

            atom_t const& atm1=parent->_mol->atom(aid1);
            int anum1=atm1.atomic_number;
            ChemData const& adata = DataForElement(anum1);
            int atomvalence=adata.nValence;
            int atomoct=adata.maxOct;

            atomoct*=0.5;
            // There are expanded octets for this system
            if(atomoct>4){
                /* only allow hypervalent atom if it has more than 2 bonds */
                if (nbonds<=2){
                    atomoct=4; 
              //}else{
                    //_component_has_expanded_octets=true;
                }
            }
            rowdata.at(apair.second)=1;
        
            BOOST_FOREACH(Id bid, bonds){
                rowdata.at(asserted_find(_component_bond_cols,bid))=1;
            }

            if(nbonds==0){
                /* Do Nothing (Dont add octet or charge constraints for ions) */
            }else if(atomoct>4){
                /* Fixes for 3-bond sulfur (probably needed for all hypervalents)...
                   only allow 10 electrons around 3 bond sulfur */
                if(anum1==16 && nbonds == 3){
                    atomoct-=1; 
                    std::vector<double> rowcopy(rowdata);
                    rowcopy.at(apair.second)=-1;
                    lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_GE,2);
                } 
                /* bounded constraint if atom can have expanded octets */
                lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_GE,4);
                lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_LE,atomoct);
            }else{
                std::vector<double> rowcopy(rowdata);
                ilpMap::iterator hiter=_component_atom_hextet_cols.find(aid1);
                if(hiter != _component_atom_hextet_cols.end()){
                    rowcopy.at(hiter->second)=1;
                }
                /* equality constraint */
                lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_EQ,atomoct);
            }
            
            ilpMap::iterator qiter=_component_atom_charge_cols.find(aid1);
            if ( qiter != _component_atom_charge_cols.end() ){
                rowdata.at(apair.second)+=1;
                /* Negative charge constraint */
                int qcol=qiter->second;
                rowdata.at(qcol)=-1;
                rowdata.at(qcol+1)=0;
                lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_LE,atomvalence);
            
                /* positive charge constraint */
                for(int k=0;k<ncols;++k) rowdata.at(k+1)*=-1;
                rowdata.at(qcol)=0; 
                rowdata.at(qcol+1)=-1;   
                lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_LE,-atomvalence);
            }
        }
    }


    void ComponentAssigner::add_component_electron_constraint(){
        //static desres::profiler::Symbol _("ComponentAssigner::add_component_electron_constraint");
        //desres::profiler::Clock __(_);

        BondOrderAssignerPtr parent=_parent.lock();
        int ncols=lpsolve::get_Norig_columns(_component_lp);
        std::vector<double> rowdata(ncols+1,0);
            
        int valence=0;
        int extvalence=0;
        BOOST_FOREACH(ilpMap::value_type const& apair, _component_atom_cols){
            Id aid = apair.first;
            int acol=apair.second;
            
            atom_t const& atm=parent->_mol->atom(aid);
            int anum=atm.atomic_number;
            ChemData const& adata = DataForElement(anum);
            valence+=adata.nValence;
            
            rowdata.at(acol)=2;
            
            IdList bonds=filteredBondsForAtom(parent->_mol, aid);
            BOOST_FOREACH(Id bid, bonds){
                rowdata.at(asserted_find(_component_bond_cols,bid))+=1;
                Id other=parent->_mol->bond(bid).other(aid);
                if (_component_atoms_present.count(other)) continue;
                /* This bonded atom is not in the component... we need to
                   adjust the component valence to take this into account */
                electronRange const& range= asserted_find(parent->_bond_order, bid);
                assert(range.lb==range.ub);
                extvalence+=range.lb;
            }
        }
        
        _component_valence_count=valence+extvalence;
     
        /* Negative charge constraint */
        rowdata.at(_component_charge_col)=-1;
        rowdata.at(_component_charge_col+1)=0;
        lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_LE,valence);
         _component_charge_row=lpsolve::get_Norig_rows(_component_lp);

        /* positive charge constraint */
        for(int k=0;k<ncols;++k) rowdata.at(k+1)*=-1;
        rowdata.at(_component_charge_col)=0; 
        rowdata.at(_component_charge_col+1)=-1;   
        lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_LE,-valence);
    }


    void ComponentAssigner::build_integer_linear_program(){
        
        //static desres::profiler::Symbol _("ComponentAssigner::build_integer_linear_program");
        //desres::profiler::Clock __(_);

        lpsolve::delete_lp(_component_lp);
        lpsolve::delete_lp(_component_lpcopy);
        _component_lp = lpsolve::make_lp(0,0);
        lpsolve::set_scaling(_component_lp,SCALE_NONE);
        /* NOTE: *ALWAYS* use at least PRESOLVE_COLS. It allows for simple determination of
           possibly resonant systems, and significantly reduces the model size. */
        int presolvetype=PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP | PRESOLVE_MERGEROWS;
        presolvetype |= PRESOLVE_IMPLIEDSLK | PRESOLVE_REDUCEGCD | PRESOLVE_BOUNDS;
        lpsolve::set_presolve(_component_lp,presolvetype, lpsolve::get_presolveloops(_component_lp) );

#if DEBUGPRINT==0
        lpsolve::set_verbose(_component_lp,0);
#endif

        // Step 1) Add columns and set objective function to minimize
        set_atom_lonepair_penalties();
        set_bond_penalties();
        set_atom_charge_penalties();
        set_atom_hextet_penalties();
        set_aromatic_ring_penalties();
        set_component_charge_penalty();

        /* The objective values for presolved columns get removed, 
           so we copy the total objective vector and save it here */
        _component_objf.assign(lpsolve::get_Norig_columns(_component_lp)+1,0);
        lpsolve::get_row(_component_lp,0,&_component_objf[0]);
        _component_objf[0]=0.0;

        /* This purturbs the objective function for the initial solution,
         but dosent affect the determination of resonance structures since we
         saved the pristine objective above */
        break_ring_symmetry();

        // Step 2) Add rows (constraints/equalities)
        lpsolve::set_add_rowmode(_component_lp, true );

        add_atom_octet_and_charge_constraints();
        add_aromatic_ring_constraints();
        add_component_electron_constraint();

        lpsolve::set_add_rowmode(_component_lp, false);

        /* after calling lpsolve::solve when presolve is in effect, its impossible to change
           any of the constraints (eg total charge) and re-solve. Therefore, we copy the model
           and update the constraints before (re)solving the model */
        _component_lpcopy=lpsolve::copy_lp(_component_lp);
    }


    /* reset status of BondOrderAssigner to initial state */
    void BondOrderAssigner::reset(){
        _valid=false;
        _total_charge_set=false;
        _total_charge=0;
        for(Id component=0; component< _component_assigners.size(); ++component){
            _component_assigners[component]->reset();
        }
    }

    void BondOrderAssigner::setTotalCharge(int qTotal){
        if((_totalValence+qTotal)%2 != 0){
            throw std::runtime_error("Desired charge results in unpaired electrons."
                                     " Only closed shell configurations are supported");
        }
        _valid=false;
        _total_charge_set=true;
        _total_charge=qTotal;
    }

    void BondOrderAssigner::unsetTotalCharge(){
        _valid=false;
        _total_charge_set=false;
        _total_charge=0;
    }
    int BondOrderAssigner::getSolvedTotalCharge(){
        if(!_valid){
            std::stringstream msg;
            msg << "Cannot getSolvedTotalCharge from system"
                " with invalid integer linear program solution." <<
                " Did you call solveIntegerLinearProgram first?";
            throw std::runtime_error(msg.str());
        }
        int qsum=_presolved_charge;
        for(Id component=0; component< _component_assigners.size(); ++component){
            qsum+=getSolvedComponentCharge(component);
        }
        if(_total_charge_set) assert(qsum==_total_charge);
        return qsum;
    }

    void BondOrderAssigner::setComponentCharge(Id component, int qTotal){
        assert(component<_component_assigners.size());
        _component_assigners[component]->setComponentCharge(qTotal);
        _fixed_component_charges[component]=qTotal;
    }
    void BondOrderAssigner::unsetComponentCharge(Id component){
        assert(component<_component_assigners.size());
        _component_assigners[component]->unsetComponentCharge();
        _fixed_component_charges.erase(component);
    }
    int BondOrderAssigner::getSolvedComponentCharge(Id component){
        if(!_valid){
            std::stringstream msg;
            msg << "Cannot getSolvedComponentCharge from component "<< component <<
                " with invalid integer linear program solution." <<
                " Did you call solveIntegerLinearProgram first?";
            throw std::runtime_error(msg.str());
        }
        assert(component<_component_assigners.size());
        int q1=_component_assigners[component]->getSolvedComponentCharge();
        std::map<Id, int>::const_iterator it = _fixed_component_charges.find(component);
        if(it != _fixed_component_charges.end()) assert(q1==it->second);
        return q1;
    }



    bool BondOrderAssigner::solveIntegerLinearProgram(){
        //static desres::profiler::Symbol _(" BondOrderAssigner::solveIntegerLinearProgram");
        //desres::profiler::Clock __(_);

        if(_needRebuild) rebuild();

        std::set<Id> active;
        int qfixed=_presolved_charge;
        for(Id component=0; component< _component_assigners.size(); ++component){
            std::map<Id, int>::const_iterator it = _fixed_component_charges.find(component);
            if(it != _fixed_component_charges.end()){
                qfixed+=it->second;
            }else{
                active.insert(component);
            }
        }

        _valid=true;
        if(_total_charge_set){
            if( active.size() == 0 ){
                _valid=(_total_charge==qfixed);
            }else if(active.size()==1){
                Id cid=*active.begin();
                _component_assigners[cid]->setComponentCharge(_total_charge-qfixed);
            }
        }
        /* try and find an initial solution */
        BOOST_FOREACH(ComponentAssignerPtr ca, _component_assigners){
            _valid &= ca->solveComponentIntegerLinearProgram();
        }

        /* Three of 4 cases are determined
         * 1) Invalid solution 
         * 2) valid solution && !_total_charge_set
         * 3) valid solution && _total_charge_set && active.size() <2
         */
        if (!_valid || !_total_charge_set || (_total_charge_set && active.size() < 2)) return _valid;
    
        /* If we get here then (valid solution && _total_charge_set && active.size()>=2).
         * Now we need to generate alternative solutions for the active components for use in
         * the total charge solver 
         */
        
        typedef Id                                      key_type;
        typedef std::pair<int,double>                   entry_type;
        typedef std::vector<entry_type>                 mapped_type;
        typedef std::pair<const key_type, mapped_type>  value_type;
        typedef std::map<key_type,mapped_type>          LpInfo;

        LpInfo solutions;

        /* Find reasonable # of valid component solutions */
        static const size_t _ndelta=4;
        static const int deltas[_ndelta]={-2,2,-4,4};
        size_t ncols=0;
        BOOST_FOREACH(Id cid, active){
            int q0=_component_assigners[cid]->getSolvedComponentCharge();
            double obj=_component_assigners[cid]->getSolvedComponentObjective();
            solutions[cid].push_back(entry_type(q0,obj));
            ncols++;
            for (size_t i=0;i<_ndelta;++i){
                int qtarget=q0+deltas[i];
                _component_assigners[cid]->setComponentCharge(qtarget);
                if(!_component_assigners[cid]->solveComponentIntegerLinearProgram()) continue;
                obj=_component_assigners[cid]->getSolvedComponentObjective();
                solutions[cid].push_back(entry_type(qtarget,obj));
                ncols++;
            }
        }

        /* Create total charge lp model */
        lpsolve::_lprec *qtotlp = lpsolve::make_lp(0,0);
        lpsolve::set_scaling(qtotlp,SCALE_NONE);
        int presolvetype=PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP | PRESOLVE_MERGEROWS;
        presolvetype |= PRESOLVE_IMPLIEDSLK | PRESOLVE_REDUCEGCD | PRESOLVE_BOUNDS;
        lpsolve::set_presolve(qtotlp, presolvetype, lpsolve::get_presolveloops(qtotlp) );
#if DEBUGPRINT==0
        lpsolve::set_verbose(qtotlp,0);
#endif

        std::ostringstream ss;
        std::vector<std::vector<double> > cons;
        std::vector<double> globalCons(ncols+1,0);
        BOOST_FOREACH(value_type const& cadata, solutions){
            cons.push_back(std::vector<double>());
            std::vector<double> &rowdata=cons.back();
            rowdata.resize(ncols+1,0);
            Id caidx=cadata.first;
            /* add new columns */
            BOOST_FOREACH(entry_type const& soldata, cadata.second){
                ss.str("");
                ss << "ca_"<<caidx<<".q_"<<soldata.first;
                int cid=ComponentAssigner::add_column_to_ilp(qtotlp,ss.str(), soldata.second, 0, 1);
                lpsolve::set_binary(qtotlp, cid, true); // redundant?
                globalCons.at(cid)=soldata.first;
                rowdata.at(cid)=1;
            }
        }

        lpsolve::set_add_rowmode(qtotlp, true );
        lpsolve::add_constraint(qtotlp, &globalCons[0], ROWTYPE_EQ, _total_charge - qfixed);
        BOOST_FOREACH(std::vector<double> &rowdata, cons){
            lpsolve::add_constraint(qtotlp, &rowdata[0], ROWTYPE_EQ, 1);
        }
        lpsolve::set_add_rowmode(qtotlp, false);
        int status=lpsolve::solve(qtotlp);
        _valid=(status==OPTIMAL || status==PRESOLVED);

        if(_valid){
#if DEBUGPRINT
            printf("Solution was found for Total Charge= %d\n", _total_charge);
            lpsolve::write_LP(qtotlp,stdout);
            lpsolve::print_objective(qtotlp);               
            lpsolve::print_solution(qtotlp,1);
#endif
            std::vector<int> solution;
            ComponentAssigner::get_ilp_solution(qtotlp,solution);
            unsigned cid=1;
            BOOST_FOREACH(value_type const& cadata, solutions){
                Id caidx=cadata.first;
                int nset=0;
                BOOST_FOREACH(entry_type const& soldata, cadata.second){
                    int sol=solution.at(cid);
                    if(sol!=0){
                        assert(sol==1);
                        _component_assigners[caidx]->setComponentCharge(soldata.first);
                        nset++;
                    }
                    cid++;
                }
                assert(nset==1);
                assert(_component_assigners[caidx]->solveComponentIntegerLinearProgram());
            }
        }

        lpsolve::delete_lp(qtotlp);
        return _valid;
    }

    double BondOrderAssigner::getSolvedObjective(){
        if(!_valid){
            throw std::runtime_error("Cannot getSolvedObjective from invalid integer linear program solution."
                                     " Did you call solveIntegerLinearProgram first?");
        }
        double obj=0.0;
        BOOST_FOREACH(ComponentAssignerPtr ca, _component_assigners){
            obj+= ca->getSolvedComponentObjective();
        }
        return obj;
    }

    void BondOrderAssigner::assignSolutionToAtoms(){
        //static desres::profiler::Symbol _("BondOrderAssigner::assignSolutionToAtoms"); 
        //desres::profiler::Clock __(_);

        if(!_valid){
            throw std::runtime_error("Cannot assignSolutionToAtoms with invalid integer linear program solution."
                                     " Did you call solveIntegerLinearProgram first?");
        }
        
        solutionMap atominfo(_atominfo);
        solutionMap bondinfo(_bondinfo);
        solutionMap chargeinfo(_chargeinfo);
  
        /* Update presolved atominfo/bondinfo/chargeinfo with solved values from componentAssigners */
        BOOST_FOREACH(ComponentAssignerPtr ca, _component_assigners){
            ca->extractComponentSolution(atominfo,bondinfo,chargeinfo);
        }

        /* Did we get everyone? */
        assert(atominfo.size()==_atom_lp.size() &&
               bondinfo.size()==_bond_order.size() &&
               chargeinfo.size()==atominfo.size());
           
        /* Assign the final charges and bond orders (both formal and resonant)
           Formal Charges are given by: fc[i]= ValenceElectrons[i] - freeElectrons[i] - 0.5*Sum_j ( BondElectrons[ij] )
        */
        Id lp_col=_mol->addAtomProp("lonepair_count",  IntType);
        Id reslp_col=_mol->addAtomProp("resonant_lonepair_count", FloatType);
        BOOST_FOREACH(Id aid1, _fragatoms){
            _mol->atom(aid1).formal_charge=0;
            _mol->atom(aid1).resonant_charge=0.0;
            _mol->atomPropValue(aid1, lp_col)=0;
            _mol->atomPropValue(aid1, reslp_col)=0.0;
        }


        /* Assign bond orders and bond electron charge part here "- 0.5*Sum_j ( BondElectrons[ij] )" */
        BOOST_FOREACH(solutionMap::value_type const& bpair, bondinfo){
            Id bid=bpair.first;
            solutionValues const& bdata=bpair.second;
    
            bond_t & bond=_mol->bond(bid);
            Id aid1=bond.i;
            atom_t& atm1=_mol->atom(aid1);
            Id aid2=bond.j;
            atom_t& atm2=_mol->atom(aid2);

            /* Bond Orders */
            int order=bdata.nonresonant;
            bond.order=order;
            double resorder=bdata.resonant;
            bond.resonant_order=resorder;

            /* Charges */
            atm1.formal_charge-=order;
            atm2.formal_charge-=order;
            atm1.resonant_charge -= resorder;
            atm2.resonant_charge -= resorder;
        }

        /* Take care of the "ValenceElectrons[i] - freeElectrons[i]" part of charge here */
        BOOST_FOREACH(solutionMap::value_type const& apair, atominfo){
            Id aid=apair.first;
            solutionValues const& adata=apair.second;
            
            _mol->atomPropValue(aid, lp_col)=adata.nonresonant;
            _mol->atomPropValue(aid, reslp_col)=adata.resonant;

            atom_t& atm=_mol->atom(aid);
            int nValence=DataForElement(atm.atomic_number).nValence;

            atm.formal_charge+= nValence - 2*adata.nonresonant;

            double resq=atm.resonant_charge;
            resq += nValence - 2*adata.resonant;
            if(fabs(resq)<1E-5){
                atm.resonant_charge=0.0;
            }else{
                atm.resonant_charge=resq;
            }

            /* Assert that the model calculated charges agree with the above.
               FIXME: just replace the above calculation of charges with whats tabulated in chargeinfo
            */
            int iqtarget=0;
            double dqtarget=0.0;
            solutionMap::const_iterator iter=chargeinfo.find(aid);
            if(iter != chargeinfo.end()){
                solutionValues const& qdata=iter->second;
                iqtarget=qdata.nonresonant;
                dqtarget=qdata.resonant;
            }
#if DEBUGPRINT
            printf("Charge for %u: %d %d  ( %f %f )\n",aid, atm.formal_charge, iqtarget, resq, dqtarget );
#endif
            assert(atm.formal_charge == iqtarget);
            assert(fabs(resq - dqtarget)<1E-8);
            
        }
    }

    int double_to_int(double dval, double tol){
        int ival=static_cast<int>(dval+(dval<0 ? -0.5 : 0.5));
        double delta=fabs(static_cast<double>(ival)-dval);
        double abstol=fabs(tol);
        if(delta>abstol){
            printf("  WARNING: Inexact integer conversion diff = %e   tol = %e\n",
                   delta,abstol);
        }
        return ival;
    }


}}
