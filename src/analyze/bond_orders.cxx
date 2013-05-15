#include <algorithm>  // sort
#include <cmath>     // fabs, pow, round
#include <stdexcept>
#include <sstream>
#include <cstdio>
#include <deque>

#include "bond_orders.hxx"
#include "../sssr.hxx"
#include "aromatic.hxx"
#include "get_fragments.hxx"
#include "../quadsum.hxx"
#include "eigensystem.hxx"
#include "../elements.hxx"

using namespace desres::msys;

namespace lpsolve {
    /* for solving the integer linear equations. We put this in a new 
       namespace since it dosent already use one */ 
#include "lp_solve/lp_lib.h"
}


namespace {
    double clamp(double v){
        static const double increment=pow(2,20); // 20 binary digits = 1048576
        return nearbyint(v*increment)/increment;
    }

    /* "Natural" electronegativity for lp */
    static const double enegLP=clamp(0.56);

    /* simple helper function to extract a value from a map (and assert that its found) */
    template<typename C>
    typename C::mapped_type const& asserted_find(C const& container, typename C::key_type const& key){
        typename C::const_iterator iter=container.find(key);
        if(iter == container.end()){
            std::stringstream msg;
            msg << "BondOrderAssigner: Programming error in call to asserted_find. " 
                << "Specified key was not found in container: " << key;
            throw std::runtime_error(msg.str());
        }
        return iter->second;
    }

    template<typename C>
    void print_Map(C const& m, std::string const& s){
        BOOST_FOREACH(typename C::value_type const& p, m){
            std::cout << s.c_str() << ": " << p.first << " " << p.second << std::endl;
        }
    }


    static inline int double_to_int(double dval) {
        int ival=static_cast<int>(dval+(dval<0 ? -0.5 : 0.5));
#if DEBUGPRINT
        double delta=fabs(static_cast<double>(ival)-dval);
        if(delta>0.0){
            printf("  WARNING: Inexact integer conversion diff = %e   tol = 0.0\n",
                   delta);
        }
#endif
        return ival;
    }

    /* Helper function for adding a column to the integer linear program
       Sets bounds and objective function contribution (as well as column 
       label for debugging) */
    int add_column_to_ilp(lpsolve::_lprec *lp, std::string const& colname, 
                                             double penalty, double lb, double ub){
        double coldata[1]={penalty};
        lpsolve::add_column(lp,coldata);
        int colid=lpsolve::get_Norig_columns(lp);
        assert(colid>0);
        if(lb==0 && ub==1){
            lpsolve::set_binary(lp,colid,true);
        }else{
            lpsolve::set_bounds(lp,colid,lb,ub);
            lpsolve::set_int(lp,colid,true);
        }
        lpsolve::set_col_name(lp,colid,const_cast<char*>(colname.c_str()));
        return colid;
    }

    void get_ilp_solution(lpsolve::_lprec *lp, std::vector<int> &solution){

        int nrows=lpsolve::get_Norig_rows(lp);
        int ncols=lpsolve::get_Norig_columns(lp);
        unsigned solsize=ncols+1;

        solution.assign(solsize,0);
        for (unsigned idx=1; idx<solsize;++idx){
            double result=lpsolve::get_var_primalresult(lp, nrows + idx);
            solution[idx]=double_to_int(result);
        }
    }

    double get_ilp_objective(std::vector<double> const& objf, 
                             std::vector<int> const& solution){
        assert(solution.size()==objf.size());
        Quadsum obj=0.0;
        for (unsigned idx=1; idx<solution.size();++idx){
            obj += objf[idx]*solution[idx];
        }
        return obj.result();
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
        boa->_filter=new bondedVirtualsAndMetalsFilter(sys);
        
        BOOST_FOREACH(Id aid, fragment){
            assert(boa->_mol->hasAtom(aid));
            boa->_mol->atom(aid).formal_charge=0; 
            boa->_mol->atom(aid).resonant_charge=0.0;
            BOOST_FOREACH(Id bid, boa->_mol->bondsForAtom(aid)){
                msys::bond_t &bond=boa->_mol->bond(bid);
                /* Initially, set kept bond orders=1.0, removed=0.0 */
                if((*boa->_filter)(bond)){
                    bond.order=1;
                    bond.resonant_order=1.0;
                }else{
                    bond.order=0;
                    bond.resonant_order=0.0;                    
                }
            }
            if(boa->_mol->atom(aid).atomic_number<1)
                //printf("BondOrderAssigner::create - Skipping atomid=%u with atomic_number<1\n",aid);
                continue;
            boa->_fragatoms.push_back(aid);
        }
        
        MultiIdList rings = GetSSSR(boa->_mol, boa->_fragatoms, true);
        /* Keep only rings with <8 atoms and all atoms <4 bonds */
        MultiIdList keepRings;
        std::vector<std::set<Id> > keepBonds;
        BOOST_FOREACH(IdList & ring, rings){
            if(ring.size()>7)continue;
            bool good=true;
            BOOST_FOREACH(Id aid, ring){
                if(boa->_mol->filteredBondedAtoms(aid, *boa->_filter).size()>3){
                    good=false;
                    break;
                }
            }
            if(good){
                keepRings.push_back(ring);
                boa->_rings.push_back(IdList());
                IdList &rbonds=boa->_rings.back();
                for (unsigned i = 0; i < ring.size(); ++i){
                    boa->_ringAtoms.insert(ring[i]);
                    Id bond = boa->_mol->findBond(ring[i],ring[(i+1)%ring.size()]);
                    if (bond == msys::BadId) MSYS_FAIL("Ring bond not found in system");
                    rbonds.push_back(bond);
                }
            }
        }
        /* this is now a MultiIdList of [rings[keepRingIds]] */
        rings=RingSystems(boa->_mol, keepRings);
        BOOST_FOREACH(IdList const& ring, rings){
            if(ring.size()==1) continue;
            std::set<Id> touched;
            BOOST_FOREACH(Id rid, ring){
                BOOST_FOREACH(Id bid, boa->_rings[rid]){
                    touched.insert(bid);
                }
            }
            boa->_rings.push_back(IdList(touched.begin(),touched.end()));
        }    


        boa->_totalValence=0;
        electronRange tmprange;
        std::set<Id> boostLP;
        BOOST_FOREACH(Id aid0, boa->_fragatoms){
            atom_t& atm0=boa->_mol->atom(aid0);
            int anum0=atm0.atomic_number;

            ChemData const& adata = DataForElement(anum0);
            IdList bonds=boa->_mol->filteredBondsForAtom(aid0, *boa->_filter);

            if(bonds.size()!=0){
                if(adata.nodata()){
                    std::cout << "Warning: No property information available for"
                              << " Atomic Number "<< anum0 
                              <<" with "<< bonds.size()<<" bonds"<<std::endl;;
                }
                if(bonds.size()>adata.maxCoord){
                    std::stringstream msg;
                    msg << "BondOrderAssigner: Something is wrong with atom " << aid0 
                        << " - Atom has " << bonds.size() 
                        << " connections. Max allowed for element " 
                        << AbbreviationForElement(anum0)
                        << " is "<< adata.maxCoord << std::endl;
                    throw std::runtime_error(msg.str());
                }
            }

            boa->_totalValence+=adata.nValence;

            tmprange.lb=0;
            tmprange.ub=boa->max_free_pairs(aid0);
#if DEBUGPRINT1
            printf("MinMax lone pairs for Aid %u : LP Range = %d %d\n",aid0,tmprange.lb,tmprange.ub);
#endif
            boa->_atom_lp.insert(std::make_pair(aid0,tmprange));

            BOOST_FOREACH(Id bid, bonds){
                Id aid1 = boa->_mol->bond(bid).other(aid0);
                if(aid0>aid1) continue;

                int maxbo=boa->max_bond_order(aid0,aid1);

#if DEBUGPRINT1
                printf("MinMax Bond Orders for Bid %u (Aids %u %u): BO Range = %d %d \n",
                       bid,aid0,aid1,1,maxbo);
#endif

                tmprange.lb=1;
                tmprange.ub=maxbo;
                boa->_bond_order.insert(std::make_pair(bid,tmprange));

            }
        }

        BOOST_FOREACH(Id aid, boostLP){
            electronMap::iterator iter=boa->_atom_lp.find(aid);
            assert(iter!=boa->_atom_lp.end());
            if(iter->second.ub<4) iter->second.ub+=1;
        }

        IdList unsolved;
        boa->presolve_octets(unsolved);
 
        /* Fill in bondinfo with presolved values */
        BOOST_FOREACH(electronMap::value_type const& epair, boa->_bond_order){
            int val=epair.second.lb;
#if DEBUGPRINT1
            printf("After presolve bond order ranges for bond %u: %d %d\n",epair.first,val,epair.second.ub);
#endif
            if(val!=epair.second.ub) continue;
            boa->_bondinfo.insert(solutionMap::value_type(epair.first, solutionValues(val,val)));
        }

        /* Fill in atominfo/chargeinfo with presolved values */
        boa->_presolved_charge=0;
        BOOST_FOREACH(electronMap::value_type const& epair, boa->_atom_lp){
            int val=epair.second.lb;
#if DEBUGPRINT1
            printf("After presolve atom lp ranges for atom %u: %d %d\n",epair.first,val,epair.second.ub);
#endif
            if(val!=epair.second.ub) continue;
            Id aid=epair.first;
            boa->_atominfo.insert(solutionMap::value_type(aid, solutionValues(val,val)));

            /* Check to see if atom was completly determined and compute charge if yes. 
             * Formal Charges are given by:  
             *    fc[i]= ValenceElectrons[i] - freeElectrons[i] - 0.5*Sum_j ( BondElectrons[j] )
             * 'val' is free electron PAIRS, so multiply by 2 */
            int qtot=DataForElement(boa->_mol->atom(aid).atomic_number).nValence - 2*val;
            bool good=true;
            BOOST_FOREACH(Id const& bid, boa->_mol->filteredBondsForAtom(aid, *boa->_filter)){
                solutionMap::const_iterator iter=boa->_bondinfo.find(bid);
                if(iter == boa->_bondinfo.end()){
                    good=false;
                    break;
                }else{
                    /* 'nonresonant' is the bond electron PAIRS, so no factor of 0.5 needed */
                    qtot-=iter->second.nonresonant;
                }
            }
            if(good){
                boa->_chargeinfo.insert(solutionMap::value_type(epair.first, solutionValues(qtot,qtot)));
                boa->_presolved_charge+=qtot;
            }
        }

        /* Split into components */
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

#if DEBUGPRINT
        printf("AfterPresolve: total_valence=%d for %zu atoms (%zu finalized,  q_of_finalized=%d)\n",
               boa->_totalValence, boa->_fragatoms.size(), boa->_atominfo.size(), boa->_presolved_charge);

#endif

        return boa;
   }

    void BondOrderAssigner::rebuild(){
        
        BOOST_FOREACH(ComponentAssignerPtr ca, _component_assigners){
            ca->build_integer_linear_program();
        }
        
        _needRebuild=false;
    }
    
    /* Cap LP counts for each atom based on group and number of bonds */
    int BondOrderAssigner::max_free_pairs(const Id aid){
        
        int maxFree;
        int nbonds=_mol->filteredBondsForAtom(aid,*_filter).size();
        int anum=_mol->atom(aid).atomic_number;
        int group=GroupForElement(anum);
        int period=PeriodForElement(anum);
        if(anum<3){
            /* hydrogen and helium */
            maxFree=1-nbonds;
        }else if(group>=1 && group <=12){
            /* metals / transition metals */
            maxFree=0;
        }else{
            /* Everything else */
            maxFree=4-nbonds;
            if(period>1 && nbonds>2)maxFree=std::max(maxFree+1,1);
        }
        return std::max(0,maxFree);
    }
    

    /* Cap bond orders for each bond connected to an atom */
    int BondOrderAssigner::max_bond_order(const Id aid0, const Id aid1){
        
        int maxOrders[2];
        Id ids[2]={aid0,aid1};
        for(int idx=0;idx<2;++idx){
            Id aid=ids[idx];  
            Id nbonds=_mol->filteredBondsForAtom(aid,*_filter).size();
            int anum=_mol->atom(aid).atomic_number;

            /* now determine max bond order */
            int maxbo;
            switch (anum){
            case 1:  // H
            case 9:  // F
                // single bond max to hydrogen, flourine
                maxbo=1;
                break;
                // other halogens (per[chlor,brom,iod]ates)
            case 17: // Cl
            case 35: // Br
            case 53: // I
                maxbo=2;
                break;
            case 5: // B
            case 6: // C
            case 7: // N
                if(nbonds == 4){
                    // single bond max to boron-, saturated carbon, nitrogen+
                    maxbo=1;
                }else if(nbonds == 3){
                    // double bond max to boron-, carbon, nitrogen+ w/3 bonds
                    maxbo=2;
                }else{
                    maxbo=3;
                }
                break;
            case 8:  // O
                if(nbonds == 1){
                    maxbo=3;
                }else{
                    /* allows [O+](-R)(=R)  */
                    maxbo=2;
                }
                break;
                
            default:
                /* Catch all... Should be fairly conservative */
                if(nbonds<3){
                    maxbo=3;
                }else if (nbonds<6){
                    maxbo=2;
                }else{
                    maxbo=1;         
                }
            }
            maxOrders[idx]=maxbo;
        }
        return std::min(maxOrders[0],maxOrders[1]);
       
    }

    /* Function that trys to presolve for unknown bond orders / lp counts
       Its can only presolve for cases where there is a single unknown.
    */
    void BondOrderAssigner::presolve_octets(IdList &unsolved){

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
                int period=PeriodForElement(atm1.atomic_number);
                if(period>2){
                    keep.push_back(aid1);
                    continue;
                }
                int octval=period<=1 ? 1 : 4; 
                
                IdList bonds=_mol->filteredBondsForAtom(aid1, *_filter);
                int unkcount=0;
                electronMap::iterator lastatom=_atom_lp.end();
                electronMap::iterator lastbond=_bond_order.end();
                
                electronMap::iterator aiter=_atom_lp.find(aid1);
#if DEBUGPRINT1
                printf("Presolving around atom %u: octval= %d   lonepair lb= %d  ub= %d\n",
                       aid1, octval, aiter->second.lb, aiter->second.ub);
#endif      
                if(aiter->second.lb==aiter->second.ub){
                    octval-=aiter->second.lb;
                }else{
                    unkcount++;
                    lastatom=aiter;
                } 
                
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
#if DEBUGPRINT1
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
#if DEBUGPRINT1
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
        lpsolve::delete_lp(_component_reslp);
    }

     ComponentAssignerPtr ComponentAssigner::create(BondOrderAssignerPtr b, IdList const& comp, Id cid){

        ComponentAssignerPtr ca(new ComponentAssigner);

        assert(b!=NULL && comp.size()>0 );
        ca->_parent=b;
        ca->_component_atoms_present=std::set<Id>(comp.begin(), comp.end());
#if DEBUGPRINT
        printf("Component %u has atoms:",cid);
        BOOST_FOREACH(Id aid, ca->_component_atoms_present) printf(" %u",aid);
        printf("\n");
#endif

        ca->_component_id=cid;
        ca->_component_valence_count=0;
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

        if(_component_solution_valid) return _component_solution_valid;
        
        lpsolve::delete_lp(_component_reslp);
        _component_reslp=lpsolve::copy_lp(_component_lp);

        int status=lpsolve::solve(_component_lp);
        _component_solution_valid=(status==OPTIMAL || status==PRESOLVED);

        if(_component_solution_valid){
#if DEBUGPRINT
            int qTotal=getSolvedComponentCharge();
            printf("Solution was found for component %u with Charge= %d   objf= %6.3f\n",
                   _component_id, qTotal, lpsolve::get_objective(_component_lp));
#endif
#if DEBUGPRINT2
            lpsolve::write_LP(_component_lp,stdout);
            lpsolve::print_objective(_component_lp);               
            lpsolve::print_solution(_component_lp,1);
#endif
        }else{
#if DEBUGPRINT
            printf("No solution found for component %u\n",_component_id);
#endif
#if DEBUGPRINT2
            lpsolve::write_LP(_component_lp,stdout);
#endif
        }

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



    void ComponentAssigner::generate_resonance_forms(){
        int ncols=lpsolve::get_Norig_columns(_component_lp);
        unsigned nvars=ncols+1;

        std::vector<int> last_solution;
        get_ilp_solution(_component_lp, last_solution);
        double objf=get_ilp_objective(_component_objf,last_solution);

        _component_resonant_solution.assign(nvars,0);
        std::vector<double> ones(nvars,1.0);
        std::vector< std::vector<int> > newcons; 

        while(true){
            for (unsigned i=1; i<nvars;++i){
                _component_resonant_solution[i]+=last_solution[i];
            }
            newcons.push_back(std::vector<int>());
            std::vector<int> &rowdata = newcons.back();
            BOOST_FOREACH(ilpAtomMap::value_type const& kv, _component_atom_cols ){
                ilpAtom const& iatom=kv.second;
                for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                    if(last_solution[icol]) rowdata.push_back(icol);
                }
            }
            BOOST_FOREACH(ilpBondMap::value_type const& kv, _component_bond_cols ){ 
                ilpBond const& ibond=kv.second;
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB;++i,++icol){
                    if(last_solution[icol]) rowdata.push_back(icol);
                }
            }
            /* FIXME: replace with unique_ptr + custom deleter */
            lpsolve::_lprec *resLP=lpsolve::copy_lp(_component_reslp);
            set_add_rowmode(resLP,true);
            BOOST_FOREACH(std::vector<int> &rowdata, newcons){
                int nactive=rowdata.size();
                if (!add_constraintex(resLP, nactive, ones.data(), rowdata.data(), ROWTYPE_LE, nactive-1)){
                    printf("Couldnt add new constraint\n");
                    assert(false);
                }
            }
            set_add_rowmode(resLP,false);
            int status=lpsolve::solve(resLP);
            if(status==OPTIMAL || status==PRESOLVED){
                get_ilp_solution(resLP, last_solution);
                double newObjf=get_ilp_objective(_component_objf,last_solution);
#if DEBUGPRINT
                printf("Resonance Solution was found for component %u objf= %6.3f\n",
                       _component_id, newObjf);
#endif
#if DEBUGPRINT2
                lpsolve::write_LP(resLP,stdout);
                lpsolve::print_objective(resLP);               
                lpsolve::print_solution(resLP,1);
#endif
                lpsolve::delete_lp(resLP); 
                if(objf != newObjf){
#if DEBUGPRINT
                   printf("resonant solution objective is larger by %e\n",newObjf-objf);
#endif
                   break;
                }
            }else{
#if DEBUGPRINT
                printf("No resonance Solution found for component %u\n",_component_id);
#endif
#if DEBUGPRINT2
                lpsolve::write_LP(resLP,stdout);
#endif
                lpsolve::delete_lp(resLP); 
                break;
            }
        }

        double nforms=newcons.size();
#if DEBUGPRINT
        printf("Found %d resonance forms\n",(int)nforms);
#endif
        for (unsigned i=1; i<nvars;++i){
            _component_resonant_solution[i]/=nforms;
        }

    }

    void ComponentAssigner::extractComponentSolution(solutionMap &atominfo,
                                                     solutionMap &bondinfo,
                                                     solutionMap &chargeinfo){

        assert(_component_solution_valid);

        get_ilp_solution(_component_lp,_component_solution);
#if GENRESFORMS
        generate_resonance_forms();
#else   
        _component_resonant_solution.resize(_component_solution.size());
        for(Id i=0; i<_component_solution.size();++i){
            _component_resonant_solution[i]=_component_solution[i];
        }
#endif


#if DEBUGPRINT1
        for(Id i=0; i<_component_solution.size();++i){
            int v1=_component_solution.at(i);
            double v2=_component_resonant_solution.at(i);
            printf("Component Solution: %u %d %f %s\n",i,
                   v1,v2, v1==v2 ? " " : "*");     
        }
#endif

        BOOST_FOREACH(ilpAtomMap::value_type const& kv, _component_atom_cols ){
            Id aid=kv.first;
            ilpAtom const& iatom=kv.second;
            int nonres=0;
            double res=0;
            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                nonres+=_component_solution.at(icol)*i; 
                res+=_component_resonant_solution.at(icol)*i;
            }
            solutionMap::iterator iter=atominfo.lower_bound(aid);
            if(iter==atominfo.end() || atominfo.key_comp()(aid,iter->first)){ 
                atominfo.insert(iter, solutionMap::value_type(aid, solutionValues(nonres, res))); 
            }else{                                                 
                assert(nonres==iter->second.nonresonant && res==iter->second.resonant); 
            }

            if(iatom.qCol){
                /* Total charge takes up 2 columns... one for + charge and one for - */
                nonres=_component_solution.at(iatom.qCol+1)-_component_solution.at(iatom.qCol);
                res=_component_resonant_solution.at(iatom.qCol+1)-_component_resonant_solution.at(iatom.qCol);
                iter=chargeinfo.lower_bound(aid);
                if(iter==chargeinfo.end() || chargeinfo.key_comp()(aid,iter->first)){ 
                    chargeinfo.insert(iter, solutionMap::value_type(aid, solutionValues(nonres, res))); 
                }else{                                                 
                    assert(nonres==iter->second.nonresonant && res==iter->second.resonant); 
                }
              
            }
        }
        BOOST_FOREACH(ilpBondMap::value_type const& kv, _component_bond_cols ){ 
            Id aid=kv.first;
            ilpBond const& ibond=kv.second;
            int nonres=0;
            double res=0;
            for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB;++i,++icol){
                nonres+=_component_solution.at(icol)*i; 
                res+=_component_resonant_solution.at(icol)*i;
            }
            solutionMap::iterator iter=bondinfo.lower_bound(aid);
            if(iter==bondinfo.end() || bondinfo.key_comp()(aid,iter->first)){ 
                bondinfo.insert(iter, solutionMap::value_type(aid, solutionValues(nonres, res))); 
            }else{                                                 
                assert(nonres==iter->second.nonresonant && res==iter->second.resonant); 
            }
        }
    }

    /* lone pair electronegativites */
    void ComponentAssigner::set_atom_lonepair_penalties(){

        _component_atom_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();
        std::ostringstream ss;
 
        double lps=clamp(parent->atom_lone_pair_scale);

        BOOST_FOREACH(Id aid1, _component_atoms_present){
            electronRange const& range= asserted_find(parent->_atom_lp,aid1);
            ss.str("");
            ss << "a_"<<aid1 << ":"<<range.lb;

            double objv=-lps*(enegLP+clamp(DataForElement(parent->_mol->atom(aid1).atomic_number).eneg));      

            int colid=add_column_to_ilp(_component_lp,ss.str(),range.lb*objv, 0, 1);
            _component_atom_cols.insert(ilpAtomMap::value_type(aid1,ilpAtom(colid,range.lb,range.ub))); 
            for(int i=range.lb+1; i<=range.ub;++i){
                ss.str("");
                ss << "a_"<<aid1 << ":"<<i;
                add_column_to_ilp(_component_lp,ss.str(),i*objv, 0, 1);
            }
        }
    }

    /* prefer bonds between more electronegative atom pairs */
    void ComponentAssigner::set_bond_penalties(){

        _component_bond_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();

        std::ostringstream ss;
        BOOST_FOREACH(Id aid1, _component_atoms_present){
            int anum1=parent->_mol->atom(aid1).atomic_number;
            double eneg1=clamp(DataForElement(anum1).eneg);
            IdList bonds=parent->_mol->filteredBondsForAtom(aid1, *parent->_filter);
            BOOST_FOREACH(Id bid, bonds){
                electronRange const& range= asserted_find(parent->_bond_order, bid);

                ilpBondMap::iterator iter=_component_bond_cols.lower_bound(bid);
                // Have we already added this bond? 
                if(!(iter==_component_bond_cols.end() || _component_bond_cols.key_comp()(bid,iter->first))) continue; 

                Id aid2=parent->_mol->bond(bid).other(aid1);
                ss.str("");
                ss << "b_"<<aid1<<"_"<<aid2<<":"<<range.lb;
                 
                int anum2=parent->_mol->atom(aid2).atomic_number;
                double hyper=(PeriodForElement(anum1)>2 || PeriodForElement(anum2)>2)? 1.01 : 1.0;
                double objv = -(eneg1 + clamp(DataForElement(anum2).eneg));
                double factor=clamp(hyper*range.lb*(1.0+0.01*range.lb*(1-range.lb)));
                int colid=add_column_to_ilp(_component_lp,ss.str(), factor*objv ,0,1);

                _component_bond_cols.insert(iter,ilpBondMap::value_type(bid,ilpBond(colid,range.lb,range.ub)));
                for(int i=range.lb+1; i<=range.ub;++i){
                    ss.str("");
                    ss << "b_"<<aid1<<"_"<<aid2<<":"<<i;
                    factor=clamp(hyper*i*(1.0+0.01*i*(1-i)));
                    add_column_to_ilp(_component_lp,ss.str(),factor*objv, 0, 1);
                }
            }
        }
    }


    bool BondOrderAssigner::allow_hextet_for_atom(Id aid1){
        int anum1=_mol->atom(aid1).atomic_number;
        Id nbonds=_mol->filteredBondsForAtom(aid1, *_filter).size();
        /* Allow hextets for... */
        if ( 
            (  (anum1==5 || anum1==13) && (nbonds<4)               ) // Al, B
            || (anum1==6               && (nbonds==3 || nbonds==2) ) // unsaturated carbon
            /*
            || (anum1==7               && (nbonds==2 || nbonds==1) ) // unsaturated nitrogen
            || (anum1==8               && (nbonds==1             ) ) // unsaturated oxygen
            */
            ) return true;
        return false;
    }

    /* penalize atoms for having a charge */
    void ComponentAssigner::set_atom_charge_penalties(){

        std::ostringstream ss;
        BondOrderAssignerPtr parent=_parent.lock();

        static const double PosZero=0.25*clamp(DataForElement(6).eneg)+0.75*clamp(DataForElement(7).eneg);
        static const double NegZero=clamp(DataForElement(9).eneg);
        double shift=clamp(parent->atom_lone_pair_scale)*(enegLP);
        double qpp=clamp(parent->atom_plus_charge_penalty);
        double qps=clamp(parent->atom_plus_charge_scale);
        double qmp=clamp(parent->atom_minus_charge_penalty);
        double qms=clamp(parent->atom_minus_charge_scale);
        double hyper=clamp(parent->hypervalent_penalty);

        BOOST_FOREACH(ilpAtomMap::value_type & kv, _component_atom_cols){
            Id aid1=kv.first;
            /* Only set charge penalties if they havent been unambiguously determined */
            // if (parent->_chargeinfo.find(aid1) != parent->_chargeinfo.end()) continue;
           
            int anum=parent->_mol->atom(aid1).atomic_number;
            double eneg=clamp(DataForElement(anum).eneg);
            double qPlus =qpp+qps*fabs(eneg-PosZero);
            double qMinus=qmp+qms*fabs(eneg-NegZero);

            ss.str("");
            ss << "qM_"<<aid1;
            /* Only need to keep track of the first column id */
            kv.second.qCol=add_column_to_ilp(_component_lp,ss.str(),qMinus+shift,0,parent->absmax_atom_charge);
                        
            ss.str("");
            ss << "qP_"<<aid1;
            add_column_to_ilp(_component_lp,ss.str(), qPlus ,0,parent->absmax_atom_charge);

            if (PeriodForElement(anum)>2){
                ss.str("");
                ss << "hyper_"<<aid1;
                int maxHyper=std::max(0,GroupForElement(anum)-14);
                kv.second.hyperCol=add_column_to_ilp(_component_lp,ss.str(), hyper,0,maxHyper);
            }

        }
    }

    /* penalize rings for not being aromatic */
    void ComponentAssigner::set_aromatic_ring_penalties(){
        
        _component_ring_cols.clear();
        BondOrderAssignerPtr parent=_parent.lock();

        double objv=clamp(parent->aromatic_ring_penalty);

        std::ostringstream ss;
        // Try to make rings aromatic 
        for(Id ridx=0; ridx<parent->_rings.size(); ++ridx){
            bool addRing=true;
            BOOST_FOREACH(Id bid, parent->_rings[ridx]){
                msys::bond_t &bond=parent->_mol->bond(bid);
                if( _component_atoms_present.count(bond.i)==0 || _component_atoms_present.count(bond.j)==0){
                    addRing=false;
                    break;
                }
            }
            if(!addRing){
                continue;
            }

            ss.str("");
            ss << "rn_"<<ridx; // 'n' in 4n+2 pi electrons for aromaticity
            int colid=add_column_to_ilp(_component_lp,ss.str(),0,0,1000);

            /* Only need to keep track of the first column id */
            _component_ring_cols.insert(ilpRingMap::value_type(ridx,colid));

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
        ss.str("");
        ss << "qCompm_"<< _component_id; 
        /* Only need to keep track of the first column id */
        _component_charge_col=add_column_to_ilp(_component_lp,ss.str(), 
               clamp(parent->component_minus_charge_penalty),0,parent->max_component_charge);
        
        ss.str("");
        ss << "qCompp_"<< _component_id; 
        add_column_to_ilp(_component_lp,ss.str(),
               clamp(parent->component_plus_charge_penalty),0,parent->max_component_charge);
        
    }

    void ComponentAssigner::add_indicator_constraints(){
        BondOrderAssignerPtr parent=_parent.lock();

        int ncols=lpsolve::get_Norig_columns(_component_lp);
        std::vector<double> rowdata;

        MultiIdList strained;

        BOOST_FOREACH(ilpAtomMap::value_type const& kv, _component_atom_cols){
            ilpAtom const& iatom=kv.second;
            rowdata.assign(ncols+1,0);

            /* Must choose one of the available lp counts */
            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                rowdata.at(icol)=1;
            }
            lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_EQ,1);

            if(parent->_ringAtoms.count(kv.first)){
                IdList bonds=parent->_mol->filteredBondsForAtom(kv.first,*parent->_filter);
                if(bonds.size()==2){
                    strained.push_back(bonds);
                }
            }

        }
        BOOST_FOREACH(ilpBondMap::value_type const& kv, _component_bond_cols){
            ilpBond const& ibond =kv.second;
            rowdata.assign(ncols+1,0);
            /* Must choose one of the available bond orders */
            for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                rowdata.at(icol)=1;
            }
            lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_EQ,1);
        }

        BOOST_FOREACH( IdList &bids, strained){
            rowdata.assign(ncols+1,0);
            ilpBond list[]={asserted_find(_component_bond_cols,bids[0]),asserted_find(_component_bond_cols,bids[1])};
            BOOST_FOREACH(ilpBond const& ibond, list){
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)=i;
                }
            }        
            lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_LE,3);
        }
    }


    void ComponentAssigner::add_atom_octet_and_charge_constraints(){

        BondOrderAssignerPtr parent=_parent.lock();

        int ncols=lpsolve::get_Norig_columns(_component_lp);
        std::vector<double> rowdata;

        BOOST_FOREACH(ilpAtomMap::value_type const& kv, _component_atom_cols){
            rowdata.assign(ncols+1,0);
            Id aid0 = kv.first;
            ilpAtom const& iatom=kv.second;
            IdList bonds=parent->_mol->filteredBondsForAtom(aid0,*parent->_filter);
            Id nbonds=bonds.size();

            atom_t const& atm0=parent->_mol->atom(aid0);
            int anum0=atm0.atomic_number;
            int atomvalence= DataForElement(anum0).nValence;
            int period=PeriodForElement(anum0);
            int atomoct=period<=1 ? 1 : 4;// adata.maxOct/2;
  
            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                rowdata.at(icol)=i;
            }
        
            BOOST_FOREACH(Id bid, bonds){
                ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)=i;
                }
            }

            if(nbonds==0){
                /* Do Nothing (Dont add octet or charge constraints for ions) */
            }else{
                std::vector<double> rowcopy(rowdata);
                if(parent->allow_hextet_for_atom(aid0)){
                    /* bound constraint */
                    lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_LE,atomoct); 
                    lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_GE,atomoct-1);
                }else{
                    if(iatom.hyperCol){
                        rowcopy[iatom.hyperCol]=-1;
                    }
                    /* equality constraint */
                    lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_EQ,atomoct);               
                }
            }
            
            if ( iatom.qCol ){
                for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                    rowdata.at(icol)+=i;
                }
                std::vector<double> rowcopy(rowdata);
                /* Negative charge constraint */
                rowcopy.at(iatom.qCol)=-1;
                rowcopy.at(iatom.qCol+1)=0;
                lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_LE,atomvalence);
            
                /* positive charge constraint */
                for(int k=0;k<ncols;++k) rowcopy.at(k+1)*=-1;
                rowcopy.at(iatom.qCol)=0; 
                rowcopy.at(iatom.qCol+1)=-1;   
                lpsolve::add_constraint(_component_lp,&rowcopy[0],ROWTYPE_LE,-atomvalence);

            }

            /* Additional fixups, force 2 connected carbon/nitrogen adjacent to terminal carbon/nitrogen to be sp hybridized
               sum of bo=4 (either =X= or -X#)
            */
            if((anum0==6 || anum0==7) && nbonds==2){
                bool sp=false;
                BOOST_FOREACH(Id bid, bonds){
                    Id aid1=parent->_mol->bond(bid).other(aid0);
                    int anum1=parent->_mol->atom(aid1).atomic_number;
                    if((anum1!=6 && anum1!=7) || 1!=parent->_mol->filteredBondsForAtom(aid1,*parent->_filter).size()) continue;
                    sp=true;
                    break;
                }
                if(sp){
                    rowdata.assign(ncols+1,0);
                    BOOST_FOREACH(Id bid, bonds){
                        ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                        for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                            rowdata.at(icol)=i;
                        }
                    }
                    lpsolve::add_constraint(_component_lp,&rowdata[0],ROWTYPE_EQ,4);
                }
            }
        }
    }


    void ComponentAssigner::add_component_electron_constraint(){

        BondOrderAssignerPtr parent=_parent.lock();
        int ncols=lpsolve::get_Norig_columns(_component_lp);
        std::vector<double> rowdata(ncols+1,0);
            
        int valence=0;
        int extvalence=0;
        BOOST_FOREACH(ilpAtomMap::value_type const& kv, _component_atom_cols){
            Id aid = kv.first;
            ilpAtom const& iatom=kv.second;

            int anum=parent->_mol->atom(aid).atomic_number;
            ChemData const& adata = DataForElement(anum);
            valence+=adata.nValence;
            
            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                rowdata.at(icol)=2*i;
            }
            
            IdList bonds=parent->_mol->filteredBondsForAtom(aid, *parent->_filter);
            BOOST_FOREACH(Id bid, bonds){
                ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)+=i;
                }
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


    void ComponentAssigner::add_aromatic_ring_constraints(){

        BondOrderAssignerPtr parent=_parent.lock();

        int ncols=lpsolve::get_Norig_columns(_component_lp);
        std::vector<double> rowdata;

        /* colid   : n in 4n+2 pi electrons
           colid+1 : electrons need to be removed to make ring aromatic
           colid+2 : electrons need to be added   to make ring aromatic
        */
        BOOST_FOREACH(ilpRingMap::value_type const& kv, _component_ring_cols){
            Id ridx=kv.first;
            int colid =kv.second;
            
            rowdata.assign(ncols+1,0);
            // pi electrons, (4*n + 2) electrons
            double target=2.0;
            rowdata.at(colid)=-4;
            
            std::set<Id> ratoms;
            BOOST_FOREACH( Id bid, parent->_rings.at(ridx)){
                msys::bond_t &bond=parent->_mol->bond(bid);
                ratoms.insert(bond.i);
                ratoms.insert(bond.j);
                ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                /* two electrons per bond order in ring... */
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)=i*2;
                }
                /* corrected for how many are available in the pi system */
                target+=2;
            }

            /* if aromaticity can be satified using alternating double bonds, do that instead
               of including lone pairs */
            if((ratoms.size()-2)%4 !=0){
                BOOST_FOREACH( Id aid, ratoms){
                    ilpAtom const& iatom = asserted_find(_component_atom_cols, aid);
                    /* two electrons per lone pair... */
                    for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                        rowdata.at(icol)=i*2;
                    }
                    /* corrected for how many are available in the pi system */
                    target+=2*(3 - parent->_mol->filteredBondsForAtom(aid, *parent->_filter).size());
                }
            }
            
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

    void ComponentAssigner::build_integer_linear_program(){

        lpsolve::delete_lp(_component_lp);
        lpsolve::delete_lp(_component_lpcopy);
        _component_lp = lpsolve::make_lp(0,0);
        lpsolve::set_scaling(_component_lp,SCALE_NONE);
        /* NOTE: *ALWAYS* use at least PRESOLVE_COLS. It allows for simple determination of
           possibly resonant systems, and significantly reduces the model size. */
        int presolvetype=PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP | PRESOLVE_MERGEROWS;
        presolvetype |= PRESOLVE_IMPLIEDSLK | PRESOLVE_REDUCEGCD | PRESOLVE_BOUNDS;
        //presolvetype=PRESOLVE_NONE;
        lpsolve::set_presolve(_component_lp,presolvetype, lpsolve::get_presolveloops(_component_lp) );

#if DEBUGPRINT==0
        lpsolve::set_verbose(_component_lp,0);
#endif

        // Step 1) Add columns and set objective function to minimize
        set_atom_lonepair_penalties();
        set_bond_penalties();
        set_atom_charge_penalties();
        set_aromatic_ring_penalties();
        set_component_charge_penalty();

        /* The objective values for presolved columns get removed, 
           so we copy the total objective vector and save it here */
        _component_objf.assign(lpsolve::get_Norig_columns(_component_lp)+1,0);
        lpsolve::get_row(_component_lp,0,&_component_objf[0]);
        _component_objf[0]=0.0;

        // Step 2) Add rows (constraints/equalities)
        lpsolve::set_add_rowmode(_component_lp, true );

        add_indicator_constraints();
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
                int cid=add_column_to_ilp(qtotlp,ss.str(), soldata.second, 0, 1);
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
#endif
#if DEBUGPRINT2
            lpsolve::write_LP(qtotlp,stdout);
            lpsolve::print_objective(qtotlp);               
            lpsolve::print_solution(qtotlp,1);
#endif
            std::vector<int> solution;
            get_ilp_solution(qtotlp,solution);
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
        if(! ( (atominfo.size()==_atom_lp.size()) &&
               (bondinfo.size()==_bond_order.size()) &&
               (chargeinfo.size()==atominfo.size()))){
            printf("Size Mismatch!\n atominfo: %lu ?==? %lu\n bondinfo: %lu ?==? %lu\n qinfo: %lu ?==? %lu\n",
                   atominfo.size(),_atom_lp.size(),bondinfo.size(),_bond_order.size(),chargeinfo.size(),atominfo.size());
            assert(false);
        };
           
        /* Assign the final charges and bond orders (both formal and resonant)
           Formal Charges are given by: fc[i]= ValenceElectrons[i] - freeElectrons[i] - 0.5*Sum_j ( BondElectrons[ij] )
        */
        BOOST_FOREACH(Id aid1, _fragatoms){
            _mol->atom(aid1).formal_charge=0;
            _mol->atom(aid1).resonant_charge=0.0;
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

#if DEBUGPRINT
            printf("Bond Order for bond %u ( %u %u ): %d %f\n",bid, aid1, aid2, order, resorder );
#endif

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
            printf("Charge for atom %u: %d %d  ( %f %f )\n",aid, atm.formal_charge, iqtarget, resq, dqtarget );
#endif
            if((atm.formal_charge != iqtarget) || fabs(resq - dqtarget)>1E-8 ){
                throw std::runtime_error("SolutionToAtoms failed. atm.formal_charge != iqtarget or fabs(resq - dqtarget)>1E-8");
            }
        }
    }
}}
