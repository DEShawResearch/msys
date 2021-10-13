#include <algorithm>  // sort
#include <cmath>     // fabs, pow, round
#include <stdexcept>
#include <sstream>
#include <cstdio>
#include <deque>
#include <iostream>
#include <memory>
#include <unordered_map>
#include <cassert>
#include <numeric>

#include "bond_orders.hxx"
#include "../sssr.hxx"
#include "get_fragments.hxx"
#include "eigensystem.hxx"
#include "../elements.hxx"

using namespace desres::msys;

namespace lpsolve {
    /* for solving the integer linear equations. We put this in a new
       namespace since it doesn't already use one */
#if defined __has_include
    #if __has_include(<lp_solve/lp_lib.h>)
        // on some systems, the lp_solve headers are dumped right in /usr/lib/include
        // or equivalent, rather than all being inside a directory called lp_solve
        #include <lp_solve/lp_lib.h>
    #else
        #include <lp_lib.h>
    #endif
#else
    //if the compiler doesn't have __has_include, then just hope that the lpsolve
    // header is in its own directory.
    #include <lp_solve/lp_lib.h>
#endif
}


namespace {
    static const double CLAMP_INCREMENT = pow(2,20); // 20 binary digits = 1048576

    double clamp(double v) {
        return nearbyint(v*CLAMP_INCREMENT) / CLAMP_INCREMENT;
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
        for (auto const& p : m){
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
                             std::vector<int> const& solution) {
        assert(solution.size()==objf.size());
        int64_t obj = 0;
        for (unsigned idx=1; idx<solution.size(); ++idx) {
            obj += (nearbyint(objf[idx] * CLAMP_INCREMENT) * solution[idx]);
        }
        return obj / CLAMP_INCREMENT;
    }
}

namespace desres { namespace msys {
    void SolutionMap::update(Id component, Id entry, int value){
        auto result = entry_to_component.insert(std::make_pair(entry, component));
        if (not result.second) {
            // each entry can belong to only one component
            if ( result.first->second != component){
                // The presolved components can get moved to a component that needs to be solved.
                // but the stored value shouldnt change. Just check and leave it in the presolved component
                if (result.first->second == Id(-1)){
                    std::vector<int> const& vec = solutionValues[entry];
                    assert(vec.size()==1 && vec[0] == value);
                    return;
                }else{
                    std::stringstream msg;
                    msg << "SolutionMap: Something is wrong with entry "<< entry
                    <<", current component is " << component << " but previous was "<<result.first->second
                    << std::endl;
                    throw std::runtime_error(msg.str());
                }
            }
        }

        std::vector<int> &vec = solutionValues[entry];
        vec.push_back(value);
        if (vec.size() != 1){
            std::stringstream msg;
            if(!resonant_solutions){
                msg << "SolutionMap: Something is wrong with component " << component
                << ", entry " << entry<<". It has already been assigned a value"
                << " but resonant solutions are not being generated!"
                << std::endl;
                throw std::runtime_error(msg.str());
            }else if(component==Id(-1)){
                msg << "SolutionMap: Something is wrong with component " << component
                << ", entry " << entry<<". It has already been assigned a value"
                << "which shouldnt happen for presolved components"
                << std::endl;
                throw std::runtime_error(msg.str());
            }
        }
    }
    boost::optional<int> SolutionMap::find(Id entry){
        boost::optional<int> answer;
        auto result = solutionValues.find(entry);
        if (result != solutionValues.end()){
            answer = result->second[0];
        }
        return answer;
    }
    Id SolutionMap::size(){
         return entry_to_component.size();
    }

    std::unordered_map<Id, solutionResult> SolutionMap::extract(){
        std::unordered_map<Id, solutionResult> results;

        for( auto const& entry: solutionValues){
            Id entryid = entry.first;
            auto const&  vec = entry.second;
            solutionResult result;
            result.nonresonant = vec[0];
            result.resonant = std::accumulate(vec.begin(), vec.end(), 0.0)/vec.size();
            results.emplace(entryid, std::move(result));
        }
        return results;
    }

    std::unordered_map<Id, std::unordered_map<Id, std::vector<int> > > SolutionMap::generate_and_remove_resonance_groups(){
        std::unordered_map<Id, std::unordered_map<Id, std::vector<int> > > groups;
        for (auto const& kv: entry_to_component){
            groups[kv.second].emplace(kv.first, std::move(solutionValues.at(kv.first)));
        }
        entry_to_component.clear();
        solutionValues.clear();

        return groups;
    }


    BondOrderAssigner::BondOrderAssigner(SystemPtr sys,
                                         IdList const& fragment,
                                         bool compute_resonant_forms,
                                         std::chrono::milliseconds timeout)
    : _compute_resonant_charge(compute_resonant_forms),
      _atominfo(compute_resonant_forms),
      _bondinfo(compute_resonant_forms),
      _chargeinfo(compute_resonant_forms)
    {
        _needRebuild=true;
        _valid=false;
        atom_lone_pair_scale=0.95;
        atom_plus_charge_penalty=0.25;
        atom_minus_charge_penalty=0.25;
        atom_plus_charge_scale=2.00;
        atom_minus_charge_scale=1.20;
        hypervalent_penalty=2.75;
        component_plus_charge_penalty=0.25;
        component_minus_charge_penalty=0.25;
        multi_bond_scale=0.9;
        aromatic_ring_penalty=0.30;
        absmax_atom_charge=2;
        max_component_charge=6;

        BondOrderAssigner* boa = this;
        boa->_needRebuild=true;
        boa->_valid=false;
        boa->_total_charge_set=false;
        boa->_total_charge=0;
        boa->_mol=sys;
        boa->_filter=new bondedVirtualsFilter(sys);
        if (timeout > std::chrono::milliseconds(0)) {
            _deadline = std::chrono::system_clock::now() + timeout;
        } else {
            _deadline = std::chrono::time_point<std::chrono::system_clock>::max();
        }

        for (Id aid : fragment){
            assert(boa->_mol->hasAtom(aid));
            boa->_mol->atom(aid).formal_charge=0;
            for (Id bid : boa->_mol->bondsForAtom(aid)){
                msys::bond_t &bond=boa->_mol->bond(bid);
                /* Initially, set kept bond orders=1.0, removed=0.0 */
                if((*boa->_filter)(bond)){
                    bond.order=1;
                    //bond.resonant_order=1.0;
                }else{
                    bond.order=0;
                    //bond.resonant_order=0.0;
                }
            }
            if(boa->_mol->atom(aid).atomic_number<1)
                continue;
            boa->_fragatoms.push_back(aid);
        }

        MultiIdList rings = GetSSSR(boa->_mol, boa->_fragatoms, true);
        /* Keep only rings with <8 atoms and all atoms <4 bonds */
        MultiIdList keepRings;
        std::vector<std::set<Id> > keepBonds;
        for (IdList & ring : rings){
            if(ring.size()>7)continue;
            bool good=true;
            for (Id aid : ring){
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
        for (IdList const& ring : rings){
            if(ring.size()==1) continue;
            std::set<Id> touched;
            for (Id rid : ring){
                for (Id bid : boa->_rings[rid]){
                    touched.insert(bid);
                }
            }
            boa->_rings.push_back(IdList(touched.begin(),touched.end()));
        }


        boa->_totalValence=0;
        electronRange tmprange;
        std::set<Id> boostLP;
        for (Id aid0 : boa->_fragatoms){
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

            for (Id bid : bonds){
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

        for (Id aid : boostLP){
            electronMap::iterator iter=boa->_atom_lp.find(aid);
            assert(iter!=boa->_atom_lp.end());
            if(iter->second.ub<4) iter->second.ub+=1;
        }

        IdList unsolved;
        boa->presolve_octets(unsolved);

        /* Fill in bondinfo with presolved values */
        for (auto const& epair : boa->_bond_order){
            int val=epair.second.lb;
#if DEBUGPRINT1
            printf("After presolve bond order ranges for bond %u: %d %d\n",epair.first,val,epair.second.ub);
#endif
            if(val!=epair.second.ub) continue;
            boa->_bondinfo.update(-1, epair.first, val);
        }

        /* Fill in atominfo/chargeinfo with presolved values */
        boa->_presolved_charge=0;
        for (auto const& epair : boa->_atom_lp){
            int val=epair.second.lb;
#if DEBUGPRINT1
            printf("After presolve atom lp ranges for atom %u: %d %d\n",epair.first,val,epair.second.ub);
#endif
            if(val!=epair.second.ub) continue;
            Id aid=epair.first;
            boa->_atominfo.update(-1, aid, val);

            /* Check to see if atom was completly determined and compute charge if yes.
             * Formal Charges are given by:
             *    fc[i]= ValenceElectrons[i] - freeElectrons[i] - 0.5*Sum_j ( BondElectrons[j] )
             * 'val' is free electron PAIRS, so multiply by 2 */
            int qtot=DataForElement(boa->_mol->atom(aid).atomic_number).nValence - 2*val;
            bool good=true;
            for (Id bid : boa->_mol->filteredBondsForAtom(aid, *boa->_filter)){
                if(boost::optional<int> found = boa->_bondinfo.find(bid)){
                    /* 'nonresonant' is the bond electron PAIRS, so no factor of 0.5 needed */
                    qtot -= found.get();
                }else{
                    good=false;
                    break;
                }
            }
            if(good){
                boa->_chargeinfo.update(-1, epair.first, qtot);
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
            for (Id cid=0; cid<components.size(); ++cid) {
                boa->_component_assigners.push_back(
                        std::make_shared<ComponentAssigner>(
                            boa, components[cid], cid));
            }
        }

#if DEBUGPRINT
        printf("AfterPresolve: total_valence=%d for %zu atoms (%zu finalized,  q_of_finalized=%d)\n",
               boa->_totalValence, boa->_fragatoms.size(), boa->_atominfo.size(), boa->_presolved_charge);

#endif
   }

    BondOrderAssigner::~BondOrderAssigner() {
        delete _filter;
    }

    void BondOrderAssigner::rebuild(){

        for (ComponentAssignerPtr ca : _component_assigners){
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
                if(nbonds==1){
                    maxbo=1;
                }else{
                    maxbo=2;
                }
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
            for (Id aid1 : unsolved){
                atom_t const& atm1=_mol->atom(aid1);
                int period=PeriodForElement(atm1.atomic_number);
                IdList bonds=_mol->filteredBondsForAtom(aid1, *_filter);
                if(period>2 && bonds.size()>0){
                    keep.push_back(aid1);
                    continue;
                }
                int octval=period<=1 ? 1 : 4;

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

                for (Id bid : bonds){
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

    long BondOrderAssigner::seconds_until_deadline() {
        std::chrono::system_clock::time_point now = std::chrono::system_clock::now();
        long seconds = std::max(1L, std::chrono::duration_cast<std::chrono::seconds>(_deadline - now).count());
        return seconds;
    }

    ComponentAssigner::ComponentAssigner(BondOrderAssigner* b,
                                         IdList const& comp,
                                         Id cid)
        : _component_lp(nullptr, [](lpsolve::_lprec* p) {})
        , _component_lpcopy(nullptr, [](lpsolve::_lprec* p) {})
    {

        ComponentAssigner* ca = this;

        assert(b!=NULL && comp.size()>0 );
        ca->_parent=b;
        ca->_component_atoms_present=std::set<Id>(comp.begin(), comp.end());
#if DEBUGPRINT
        printf("Component %u has atoms:",cid);
        for (Id aid : ca->_component_atoms_present) printf(" %u",aid);
        printf("\n");
#endif

        ca->_component_id=cid;
        ca->_component_valence_count=0;
        ca->_component_solution_valid=false;
        ca->_component_charge_set=false;
    }

    void ComponentAssigner::reset(){
        _component_lp = unique_lp(lpsolve::copy_lp(_component_lpcopy.get()), lpsolve::delete_lp);
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
            lpsolve::set_bounds(_component_lp.get(), _component_charge_col  , -qTotal,-qTotal);
            lpsolve::set_bounds(_component_lp.get(), _component_charge_col+1, 0,      0);
            lpsolve::set_constr_type(_component_lp.get(),_component_charge_row  , ROWTYPE_EQ);
            lpsolve::set_constr_type(_component_lp.get(),_component_charge_row+1, ROWTYPE_LE);
        }else{
            lpsolve::set_bounds(_component_lp.get(), _component_charge_col  , 0,      0);
            lpsolve::set_bounds(_component_lp.get(), _component_charge_col+1, qTotal, qTotal);
            lpsolve::set_constr_type(_component_lp.get(),_component_charge_row  , ROWTYPE_LE);
            lpsolve::set_constr_type(_component_lp.get(),_component_charge_row+1, ROWTYPE_EQ);
        }
        _component_charge_set=true;
    }

    void ComponentAssigner::unsetComponentCharge(){
        if(_component_solution_valid) reset();
        BondOrderAssigner* parent=_parent;
        lpsolve::set_bounds(_component_lp.get(), _component_charge_col,  0, parent->max_component_charge);
        lpsolve::set_bounds(_component_lp.get(), _component_charge_col+1,0, parent->max_component_charge);

        lpsolve::set_constr_type(_component_lp.get(),_component_charge_row  , ROWTYPE_LE);
        lpsolve::set_constr_type(_component_lp.get(),_component_charge_row+1, ROWTYPE_LE);
        _component_charge_set=false;
    }

    int ComponentAssigner::getSolvedComponentCharge(){
        if(!_component_solution_valid){
            std::stringstream msg;
            msg << "Cannot getSolvedComponentCharge from component "<<_component_id<<
                " with invalid integer linear program solution." <<
                " Did you call solveIntegerLinearProgram first?";
            throw std::runtime_error(msg.str());
        }
        int nrows=lpsolve::get_Norig_rows(_component_lp.get());
        int qTotal=
            - double_to_int(lpsolve::get_var_primalresult(_component_lp.get(), nrows + _component_charge_col))
            + double_to_int(lpsolve::get_var_primalresult(_component_lp.get(), nrows + _component_charge_col+1));
        return qTotal;
    }

    bool ComponentAssigner::solveComponentIntegerLinearProgram(){

        if(_component_solution_valid) return _component_solution_valid;

        lpsolve::set_timeout(_component_lp.get(), _parent->seconds_until_deadline());
        int status=lpsolve::solve(_component_lp.get());
        if ((status == SUBOPTIMAL) && (std::chrono::system_clock::now() > _parent->_deadline)) {
            status = TIMEOUT;
        }
        _component_solution_valid=(status==OPTIMAL || status==PRESOLVED);

        if(_component_solution_valid){
#if DEBUGPRINT
            int qTotal=getSolvedComponentCharge();
            printf("Solution was found for component %u with Charge= %d   objf= %6.3f\n",
                   _component_id, qTotal, lpsolve::get_objective(_component_lp.get()));
#endif
#if DEBUGPRINT2
            lpsolve::write_LP(_component_lp.get(), stdout);
            lpsolve::print_objective(_component_lp.get());
            lpsolve::print_solution(_component_lp.get(), 1);
#endif
        }else{
#if DEBUGPRINT
            printf("No solution found for component %u\n",_component_id);
#endif
#if DEBUGPRINT2
            lpsolve::write_LP(_component_lp.get(), stdout);
#endif
        }

        if(!_component_solution_valid) reset();
        if (status == TIMEOUT) {
            throw std::runtime_error("Unable to solve ILP. Timeout elapsed");
        }
        if (status == INFEASIBLE) {
            throw std::runtime_error("Unable to solve ILP. Infeasible");
        }
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
        return lpsolve::get_objective(_component_lp.get());
    }

#ifdef MSYS_WITHOUT_BLISS
    std::vector<std::vector<int>> ComponentAssigner::generate_resonance_ilp_permutations(
        lpsolve::_lprec* lp,
        const std::vector<int>& solution
    ) {
        throw std::runtime_error("msys was not compiled with BLISS support");
    }
#endif

    void ComponentAssigner::generate_resonance_forms_old() {
        _resonant_solutions.clear();

        generate_resonance_forms_core([&](
            std::vector<ilp_resonance_constraint>& newcons,
            const std::vector<int>& lastSolution,
            lpsolve::_lprec* resLP
        ) {
            // http://yetanothermathprogrammingconsultant.blogspot.com/2011/10/integer-cuts.html
            ilp_resonance_constraint cons;

            for (ilpAtomMap::value_type const& kv : _component_atom_cols) {
                ilpAtom const& iatom=kv.second;
                for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                    if (lastSolution[icol]) {
                        cons.colno.push_back(icol);
                    }
                }
            }

            for (ilpBondMap::value_type const& kv : _component_bond_cols ){
                ilpBond const& ibond=kv.second;
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB;++i,++icol){
                    if (lastSolution[icol]) {
                        cons.colno.push_back(icol);
                    }
                }
            }

            cons.coeff.assign(cons.colno.size(), 1.0);
            cons.rhs = cons.colno.size() - 1.0;
            cons.rowtype = ROWTYPE_LE;

            if (cons.colno.size() == 0) {
                /* Corner Case: need to actually exclude the solution that has no set columns */
                int ncols = lpsolve::get_Norig_columns(_component_lp.get());
                unsigned nvars = ncols + 1;
                cons.colno.resize(nvars);
                cons.coeff.resize(nvars);
                cons.rhs = -1;
                for (auto i=0u; i < nvars; i++) {
                    cons.colno[i] = i+1;
                    cons.coeff[i] = 1.0;
		        }
            }

            newcons.push_back(cons);

            _resonant_solutions.push_back(lastSolution);

        });
    }

    void ComponentAssigner::generate_resonance_forms_new() {
        _resonant_solutions.clear();

        generate_resonance_forms_core([&](
            std::vector<ilp_resonance_constraint>& newcons,
            const std::vector<int>& lastSolution,
            lpsolve::_lprec* resLP
        ) {
            // generate the permutations of the last solution. We need to pass in a version
            // that has not been presolved.
            std::vector<std::vector<int>> permutations = generate_resonance_ilp_permutations(resLP, lastSolution);
            for (auto const& solution: permutations) {
                // http://yetanothermathprogrammingconsultant.blogspot.com/2011/10/integer-cuts.html
                ilp_resonance_constraint cons;
                cons.rowtype = ROWTYPE_LE;

                for (ilpAtomMap::value_type const& kv : _component_atom_cols) {
                    ilpAtom const& iatom=kv.second;
                    for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                        if (icol > 0 && solution[icol] == 1) {
                            cons.colno.push_back(icol);
                            cons.coeff.push_back(1.0);
                        }
                    }
                }

                for (ilpBondMap::value_type const& kv : _component_bond_cols ){
                    ilpBond const& ibond=kv.second;
                    for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB;++i,++icol){
                        if (icol > 0 && solution[icol] == 1) {
                            cons.colno.push_back(icol);
                            cons.coeff.push_back(1.0);
                        }
                    }
                }

                if (cons.colno.size() > 0) {
                    cons.coeff.assign(cons.colno.size(), 1.0);
                    cons.rhs = cons.colno.size() - 1.0;
                } else {
                    /* Corner Case: need to actually exclude the solution that has no set columns */
                    int ncols = lpsolve::get_Norig_columns(_component_lp.get());
                    unsigned nvars = ncols + 1;
                    cons.colno.resize(nvars);
                    cons.coeff.resize(nvars);
                    cons.rhs = -1;
                    for (auto i=0u; i < nvars; i++) {
                        cons.colno[i] = i+1;
                        cons.coeff[i] = 1.0;
                    }
                }
                newcons.push_back(cons);
            }

            assert(permutations.size() > 0);

            for (std::vector<int> &solution: permutations) {
                _resonant_solutions.push_back(std::move(solution));
            }

        });
    }

    void ComponentAssigner::generate_resonance_forms_core(
        std::function<void(
            std::vector<ilp_resonance_constraint>&,
            const std::vector<int>&,
            lpsolve::_lprec*
        )>
        build_resonance_constraints_from_soln
    ) {
        int ncols=lpsolve::get_Norig_columns(_component_lp.get());
        unsigned nvars=ncols+1;

        std::vector<int> last_solution;
        get_ilp_solution(_component_lp.get(), last_solution);
        double objf = get_ilp_objective(_component_objf,last_solution);

        std::vector<ilp_resonance_constraint> newcons;
        std::vector<ilp_resonance_constraint> accelcons;
        // these constraints are built in the original indexing, before any
        // presolves. but _component_lp might have some presolved columns, so
        // when we apply them we need to apply them to a copy of _component_lpcopy,
        // since that one was never (pre)solved
        build_acceleration_resonance_constraints(last_solution, accelcons);
        int numSolves = 1;

        while (true) {
            auto resLP = unique_lp(lpsolve::copy_lp(_component_lpcopy.get()), lpsolve::delete_lp);
            assert(1 + lpsolve::get_Norig_columns(resLP.get()) == static_cast<int>(_component_objf.size()));
            assert(lpsolve::get_Ncolumns(resLP.get()) == lpsolve::get_Norig_columns(resLP.get()));
            assert(lpsolve::get_Nrows(resLP.get()) == lpsolve::get_Norig_rows(resLP.get()));

            build_resonance_constraints_from_soln(newcons, last_solution, resLP.get());

	        assert(lpsolve::set_add_rowmode(resLP.get(), true));
            add_ilp_resonance_constraints(resLP.get(), accelcons);
            add_ilp_resonance_constraints(resLP.get(), newcons);
	        assert(lpsolve::set_add_rowmode(resLP.get(), false));

#if DEBUGPRINT2
            // better to print before solving, since presolving changes the structure
            // of the problem a bit
            lpsolve::write_LP(resLP.get(), stdout);
#endif
            lpsolve::set_timeout(resLP.get(), _parent->seconds_until_deadline());
            int status = lpsolve::solve(resLP.get());
            numSolves++;

            if ((status == SUBOPTIMAL) || (std::chrono::system_clock::now() > _parent->_deadline)) {
                status = TIMEOUT;
                throw std::runtime_error("Timeout reached");
            }
            if (status==OPTIMAL || status==PRESOLVED){
                get_ilp_solution(resLP.get(), last_solution);
                last_solution.resize(nvars);
                double newObjf = get_ilp_objective(_component_objf, last_solution);
#if DEBUGPRINT
                printf("Resonance Solution was found for component %u objf= %6.3f\n",
                    _component_id, newObjf);
#endif
#if DEBUGPRINT2
                if (status == OPTIMAL || status == PRESOLVED) {
                    lpsolve::print_objective(resLP.get());
                    lpsolve::print_solution(resLP.get(), 1);
                }
#endif
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
                break;
            }
            if (status == TIMEOUT) {
                throw std::runtime_error("Unable to solve ILP. Timeout elapsed");
            }
        }
#if DEBUGPRINT
        printf("Found %d resonance forms with %d ILP solutions\n", _resonant_solutions.size(), numSolves-1);
#endif

    }

    void ComponentAssigner::build_acceleration_resonance_constraints(
        const std::vector<int>& solution,
        std::vector<ilp_resonance_constraint>& newcons
    ) {
        // Add constraints to enforce the correct objective value
        // These are not essential for correctness, but have been observed to speed up the process
        // of finding distinct resonance solutions

        /* 1) group by variables that have the same contribution to the objective */
        std::unordered_map<int, std::vector<int> > unique_envs;
        for (unsigned i=1; i<_component_objf.size(); ++i){
            int key = nearbyint(_component_objf[i] * CLAMP_INCREMENT);
            unique_envs[key].push_back(static_cast<int>(i));
        }

        /* 2) add a constraint that enforces the same number of active variables within
           each group of equivalent objective contributions. There is no way to get the same
           objective if this is violated */
        for (auto & kv: unique_envs){
            std::vector<int> colno=kv.second;
            if (colno.size() < 2) {
                continue;
            }

            int sum = 0;
            for (int icol: colno){
                sum += solution[icol];
            }

            std::vector<double> coeff(colno.size(), 1.0);
            newcons.push_back(ilp_resonance_constraint{
                colno,
                coeff,
                static_cast<double>(sum),
                ROWTYPE_EQ
            });
        }
    }

    void ComponentAssigner::add_ilp_resonance_constraints(
        lpsolve::_lprec* lp,
        const std::vector<ilp_resonance_constraint>& constraints)
    {
        assert (lpsolve::get_Norig_columns(lp) == lpsolve::get_Ncolumns(lp));
        for (auto i = 0u; i < constraints.size(); i++) {
            auto const& c = constraints[i];
            assert(c.colno.size() == c.coeff.size());
            for (auto colid : c.colno) {
                if (!(colid > 0 && colid <= lpsolve::get_Norig_columns(lp)+1)) {
                    printf("bad colno=%d %d\n", colid, lpsolve::get_Norig_columns(lp));
                    assert(false);
                }
            }
            if (!lpsolve::add_constraintex(lp, c.colno.size(), const_cast<double*>(&c.coeff[0]), const_cast<int*>(&c.colno[0]), c.rowtype, c.rhs)) {
                fprintf(stderr, "Failed to add new constraint\n");
                assert(false);
            }
        }
    }


    void ComponentAssigner::update_infos_for_solution(std::vector<int> const& solution,
                                                         SolutionMap &atominfo,
                                                         SolutionMap &bondinfo,
                                                         SolutionMap &chargeinfo){
        for (ilpAtomMap::value_type const& kv : _component_atom_cols ){
            Id aid=kv.first;
            ilpAtom const& iatom=kv.second;
            int value=0;
            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                value+=solution.at(icol)*i;
            }
            atominfo.update(_component_id, aid, value);

            if(iatom.qCol){
                /* Total charge takes up 2 columns... one for + charge and one for - */
                value = solution.at(iatom.qCol+1) - solution.at(iatom.qCol);
                chargeinfo.update(_component_id, aid, value);
            }
        }

        for (ilpBondMap::value_type const& kv : _component_bond_cols ){
            Id aid=kv.first;
            ilpBond const& ibond=kv.second;
            int value=1;
            for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB;++i,++icol){
                value += solution.at(icol)*(i-1);
            }
            bondinfo.update(_component_id, aid, value);
        }
    }


    void ComponentAssigner::extractComponentSolution(SolutionMap &atominfo,
                                                     SolutionMap &bondinfo,
                                                     SolutionMap &chargeinfo){

        assert(_component_solution_valid);
        get_ilp_solution(_component_lp.get(), _component_solution);
        update_infos_for_solution(_component_solution, atominfo, bondinfo, chargeinfo);
        if (_parent->compute_resonant_charge()) {
            #ifdef MSYS_WITHOUT_BLISS
                generate_resonance_forms_old();
            #else
                generate_resonance_forms_new();
            #endif
            for (auto const& solution: _resonant_solutions){
                // The original _component_solution is one of the
                // resonant solutions, so we skip it
                if (solution == _component_solution){
                    continue;
                }
                update_infos_for_solution(solution, atominfo, bondinfo, chargeinfo);
            }
        }
    }

    /* lone pair electronegativites */
    void ComponentAssigner::set_atom_lonepair_penalties(){

        _component_atom_cols.clear();
        BondOrderAssigner* parent=_parent;
        std::ostringstream ss;

        double lps = parent->atom_lone_pair_scale;

        for (Id aid1 : _component_atoms_present){
            electronRange const& range= asserted_find(parent->_atom_lp,aid1);
            assert(range.lb==0);

            int anum1 = parent->_mol->atom(aid1).atomic_number;
            double eneg = enegLP + DataForElement(anum1).eneg;
            double objv = clamp(-lps*eneg);

            std::vector<int> colid;
            for(int i=1; i<=range.ub;++i){
                ss.str("");
                ss << "a_"<<aid1 << "_"<<i;
                colid.push_back(add_column_to_ilp(_component_lp.get(), ss.str(),i*objv, 0, 1));
            }
            if(colid.size())
                _component_atom_cols.insert(ilpAtomMap::value_type(aid1,ilpAtom(colid[0],1,range.ub)));
            else
                _component_atom_cols.insert(ilpAtomMap::value_type(aid1,ilpAtom(0,0,0)));

        }
    }

    /* prefer bonds between more electronegative atom pairs */
    void ComponentAssigner::set_bond_penalties(){

        _component_bond_cols.clear();
        BondOrderAssigner* parent=_parent;

        double boscale=clamp(parent->multi_bond_scale);

        std::ostringstream ss;
        for (Id aid1 : _component_atoms_present){
            int anum1=parent->_mol->atom(aid1).atomic_number;
            double eneg1=clamp(DataForElement(anum1).eneg);
            IdList bonds=parent->_mol->filteredBondsForAtom(aid1, *parent->_filter);
            for (Id bid : bonds){
                electronRange const& range= asserted_find(parent->_bond_order, bid);

                ilpBondMap::iterator iter=_component_bond_cols.lower_bound(bid);
                // Have we already added this bond?
                if(!(iter==_component_bond_cols.end() || _component_bond_cols.key_comp()(bid,iter->first))) continue;
                assert(range.lb==1);

                Id aid2=parent->_mol->bond(bid).other(aid1);

                int anum2=parent->_mol->atom(aid2).atomic_number;
                double hyper=(PeriodForElement(anum1)>2 || PeriodForElement(anum2)>2)? 1.01 : 1.0;
                double objv = -(eneg1 + clamp(DataForElement(anum2).eneg));

                std::vector<int> colid;
                double factor=1.0*hyper;
                for(int i=2; i<=range.ub;++i){
                    ss.str("");
                    ss << "b_"<<aid1<<"_"<<aid2<<"_"<<i;
                    colid.push_back(add_column_to_ilp(_component_lp.get(), ss.str(), factor*objv, 0, 1));
                    factor+=clamp(hyper*boscale);
                }
                if(colid.size())
                    _component_bond_cols.insert(iter,ilpBondMap::value_type(bid,ilpBond(colid[0],2,range.ub)));
                else
                    _component_bond_cols.insert(iter,ilpBondMap::value_type(bid,ilpBond(0,1,1)));

            }
        }
    }


    bool BondOrderAssigner::allow_hextet_for_atom(Id aid1){
        int anum1=_mol->atom(aid1).atomic_number;
        Id nbonds=_mol->filteredBondsForAtom(aid1, *_filter).size();
        /* Allow hextets for... */
        if (
            ( (anum1==5 || anum1==13) && (nbonds<4)               ) || // Al, B
            ( (anum1==6 || anum1==14) && (nbonds<4) ) || // unsaturated carbon
            ( (anum1==7 || anum1==15) && (nbonds<3) ) || // unsaturated nitrogen
            ( (anum1==8 || anum1==16) && (nbonds<2) )    // unsaturated oxygen
            ) return true;
        return false;
    }

    /* penalize atoms for having a charge (applied to all atoms in the system) */
    void ComponentAssigner::set_atom_charge_penalties(){

        std::ostringstream ss;
        BondOrderAssigner* parent=_parent;

        static const double PosZero=0.25*clamp(DataForElement(6).eneg)+0.75*clamp(DataForElement(7).eneg);
        static const double NegZero=clamp(DataForElement(9).eneg);
        double shift=clamp(parent->atom_lone_pair_scale)*(enegLP);
        double qpp=clamp(parent->atom_plus_charge_penalty);
        double qps=clamp(parent->atom_plus_charge_scale);
        double qmp=clamp(parent->atom_minus_charge_penalty);
        double qms=clamp(parent->atom_minus_charge_scale);
        double hyper=clamp(parent->hypervalent_penalty);

        for (ilpAtomMap::value_type & kv : _component_atom_cols){
            Id aid1=kv.first;

            int anum=parent->_mol->atom(aid1).atomic_number;
            double eneg=clamp(DataForElement(anum).eneg);
            double qPlus =qpp+qps*fabs(eneg-PosZero);
            double qMinus=qmp+qms*fabs(eneg-NegZero);

            ss.str("");
            ss << "qM_"<<aid1;
            /* Only need to keep track of the first column id */
            if(GroupForElement(anum)==18){
                /* nobel gases shouldnt be negative */
                kv.second.qCol=add_column_to_ilp(_component_lp.get(), ss.str(),qMinus+shift,0, 0);
            }else{
                kv.second.qCol=add_column_to_ilp(_component_lp.get(), ss.str(),qMinus+shift,0,parent->absmax_atom_charge);
            }


            ss.str("");
            ss << "qP_"<<aid1;
            if(GroupForElement(anum)>16){
                /* halogens and nobel gases shouldnt be positive */
                add_column_to_ilp(_component_lp.get(), ss.str(), qPlus ,0, 0);
            }else{
                add_column_to_ilp(_component_lp.get(), ss.str(), qPlus ,0,parent->absmax_atom_charge);
            }

            if (PeriodForElement(anum)>2){
                ss.str("");
                ss << "hyper_"<<aid1;
                int maxHyper=std::max(0,GroupForElement(anum)-14);
                kv.second.hyperCol=add_column_to_ilp(_component_lp.get(), ss.str(), hyper,0,maxHyper);
            }

        }
    }

    /* penalize rings for not being aromatic */
    void ComponentAssigner::set_aromatic_ring_penalties(){

        _component_ring_cols.clear();
        BondOrderAssigner* parent=_parent;

        double objv=clamp(parent->aromatic_ring_penalty);

        std::ostringstream ss;
        // Try to make rings aromatic
        for(Id ridx=0; ridx<parent->_rings.size(); ++ridx){
            bool addRing=true;
            for (Id bid : parent->_rings[ridx]){
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
            int colid=add_column_to_ilp(_component_lp.get(),ss.str(),0,0,1000);

            /* Only need to keep track of the first column id */
            _component_ring_cols.insert(ilpRingMap::value_type(ridx,colid));

            ss.str("");
            ss << "rs_"<<ridx; // do we need to subtract electrons from ring to make it aromatic?
            add_column_to_ilp(_component_lp.get(), ss.str(), objv, 0, 2);

            ss.str("");
            ss << "ra_"<<ridx; // or add electrons to ring to make it aromatic?
            add_column_to_ilp(_component_lp.get(), ss.str(), objv, 0, 2);

        }
    }

    /* penalize components for having a charge */
    void ComponentAssigner::set_component_charge_penalty(){

        BondOrderAssigner* parent=_parent;

        std::ostringstream ss;
        ss.str("");
        ss << "qCompm_"<< _component_id;
        /* Only need to keep track of the first column id */
        _component_charge_col=add_column_to_ilp(_component_lp.get(), ss.str(),
               clamp(parent->component_minus_charge_penalty),0,parent->max_component_charge);

        ss.str("");
        ss << "qCompp_"<< _component_id;
        add_column_to_ilp(_component_lp.get(), ss.str(),
               clamp(parent->component_plus_charge_penalty),0,parent->max_component_charge);

    }

    void ComponentAssigner::add_indicator_constraints(){
        BondOrderAssigner* parent=_parent;

        int ncols=lpsolve::get_Norig_columns(_component_lp.get());
        std::vector<double> rowdata;

        MultiIdList strained;

        for (ilpAtomMap::value_type const& kv : _component_atom_cols){
            if(parent->_ringAtoms.count(kv.first)){
                IdList bonds=parent->_mol->filteredBondsForAtom(kv.first,*parent->_filter);
                if(bonds.size()==2){
                    strained.push_back(bonds);
                }
            }
            ilpAtom const& iatom=kv.second;
            if (iatom.ilpCol==0) continue;

            rowdata.assign(ncols+1,0);
            /* Must choose at most one of the available lp counts */
            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                rowdata.at(icol)=1;
            }
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_GE,0);
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_LE,1);

        }
        for (ilpBondMap::value_type const& kv : _component_bond_cols){
            ilpBond const& ibond =kv.second;
            if (ibond.ilpCol==0) continue;
            rowdata.assign(ncols+1,0);
            /* Must choose at most one of the available bond orders */
            for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                rowdata.at(icol)=1;
            }
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_GE,0);
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_LE,1);
        }

        /* adjacent bonds cant both have bo > 1 in rings */
        for ( IdList &bids : strained){
            rowdata.assign(ncols+1,0);
            ilpBond list[]={asserted_find(_component_bond_cols,bids[0]),asserted_find(_component_bond_cols,bids[1])};
            for (ilpBond const& ibond : list){
                if (ibond.ilpCol==0) continue;
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)=(i-1);
                }
            }
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_LE,1);
        }
    }


    void ComponentAssigner::add_atom_octet_and_charge_constraints(){

        BondOrderAssigner* parent=_parent;

        int ncols=lpsolve::get_Norig_columns(_component_lp.get());
        std::vector<double> rowdata;

        for (ilpAtomMap::value_type const& kv : _component_atom_cols){
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

            for (Id bid : bonds){
                ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)=i-1;
                }
            }
            atomoct-=nbonds;
            atomvalence-=nbonds;

            if(nbonds==0){
                /* Do Nothing (Dont add octet or charge constraints for ions) */
            }else{
                std::vector<double> rowcopy(rowdata);
                if(parent->allow_hextet_for_atom(aid0)){
                    /* bound constraint */
                    lpsolve::add_constraint(_component_lp.get(), &rowcopy[0],ROWTYPE_LE,atomoct);
                    lpsolve::add_constraint(_component_lp.get(), &rowcopy[0],ROWTYPE_GE,std::max(atomoct-2,0));
                }else{
                    if(iatom.hyperCol){
                        rowcopy[iatom.hyperCol]=-1;
                    }
                    /* equality constraint */
                    lpsolve::add_constraint(_component_lp.get(), &rowcopy[0],ROWTYPE_EQ,atomoct);
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
                lpsolve::add_constraint(_component_lp.get(), &rowcopy[0],ROWTYPE_LE,atomvalence);

                /* positive charge constraint */
                for(int k=0;k<ncols;++k) rowcopy.at(k+1)*=-1;
                rowcopy.at(iatom.qCol)=0;
                rowcopy.at(iatom.qCol+1)=-1;
                lpsolve::add_constraint(_component_lp.get(), &rowcopy[0],ROWTYPE_LE,-atomvalence);

            }

            /* Additional fixups, force 2 connected carbon/nitrogen adjacent to terminal carbon/nitrogen to be sp hybridized
               sum of bo=4 (either =X=Y or -X#Y)
            */
            if((anum0==6 || anum0==7) && nbonds==2 && !_component_charge_set){
                bool sp=false;
                bool hasHydrogen=false;
                for (Id bid : bonds){
                    Id aid1=parent->_mol->bond(bid).other(aid0);
                    int anum1=parent->_mol->atom(aid1).atomic_number;
                    if((anum1==6 || anum1==7) && 1==parent->_mol->filteredBondsForAtom(aid1,*parent->_filter).size()){
                        sp=true;
                    }else if(anum1==1){
                        hasHydrogen=true;
                    }
                }
                if(sp && !hasHydrogen){
                    rowdata.assign(ncols+1,0);
                    for (Id bid : bonds){
                        ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                        for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                            rowdata.at(icol)=i-1;
                        }
                    }
                    lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_EQ,2);
                }
            }
        }
    }


    void ComponentAssigner::add_component_electron_constraint(){

        BondOrderAssigner* parent=_parent;
        int ncols=lpsolve::get_Norig_columns(_component_lp.get());
        std::vector<double> rowdata(ncols+1,0);

        int valence=0;
        int extvalence=0;
        for (ilpAtomMap::value_type const& kv : _component_atom_cols){
            Id aid = kv.first;
            ilpAtom const& iatom=kv.second;

            int anum=parent->_mol->atom(aid).atomic_number;
            ChemData const& adata = DataForElement(anum);
            valence+=adata.nValence;

            for(int icol=iatom.ilpCol, i=iatom.ilpLB; i<=iatom.ilpUB;++i,++icol){
                rowdata.at(icol)=2*i;
            }

            IdList bonds=parent->_mol->filteredBondsForAtom(aid, *parent->_filter);
            for (Id bid : bonds){
                ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)+=(i-1);
                }
                Id other=parent->_mol->bond(bid).other(aid);
                if (_component_atoms_present.count(other)) continue;
                /* This bonded atom is not in the component... we need to
                   adjust the component valence to take this into account */
                electronRange const& range= asserted_find(parent->_bond_order, bid);
                assert(range.lb==range.ub);
                extvalence+=range.lb;
            }
            valence-=bonds.size();
            extvalence+=bonds.size();
        }

        _component_valence_count=valence+extvalence;

        /* Negative charge constraint */
        rowdata.at(_component_charge_col)=-1;
        rowdata.at(_component_charge_col+1)=0;
        lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_LE,valence);
         _component_charge_row=lpsolve::get_Norig_rows(_component_lp.get());

        /* positive charge constraint */
        for(int k=0;k<ncols;++k) rowdata.at(k+1)*=-1;
        rowdata.at(_component_charge_col)=0;
        rowdata.at(_component_charge_col+1)=-1;
        lpsolve::add_constraint(_component_lp.get(), &rowdata[0],ROWTYPE_LE,-valence);
    }


    void ComponentAssigner::add_aromatic_ring_constraints(){

        BondOrderAssigner* parent=_parent;

        int ncols=lpsolve::get_Norig_columns(_component_lp.get());
        std::vector<double> rowdata;

        /* colid   : n in 4n+2 pi electrons
           colid+1 : electrons need to be removed to make ring aromatic
           colid+2 : electrons need to be added   to make ring aromatic
        */
        for (ilpRingMap::value_type const& kv : _component_ring_cols){
            Id ridx=kv.first;
            int colid =kv.second;

            rowdata.assign(ncols+1,0);
            // pi electrons, (4*n + 2) electrons
            double target=2.0;
            rowdata.at(colid)=-4;

            std::set<Id> ratoms;
            for ( Id bid : parent->_rings.at(ridx)){
                msys::bond_t &bond=parent->_mol->bond(bid);
                ratoms.insert(bond.i);
                ratoms.insert(bond.j);
                ilpBond const& ibond = asserted_find(_component_bond_cols,bid);
                /* two electrons per bond order in ring... */
                for(int icol=ibond.ilpCol, i=ibond.ilpLB; i<=ibond.ilpUB; ++i,++icol){
                    rowdata.at(icol)=(i-1)*2;
                }
            }

            /* if aromaticity can be satified using alternating double bonds, do that instead
               of including lone pairs */
            if((ratoms.size()-2)%4 !=0){
                for ( Id aid : ratoms){
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
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0], ROWTYPE_LE, target);

            /* constraint on adding electrons to ring */
            for(int k=0;k<ncols;++k) rowdata.at(k+1)*=-1;
            rowdata.at(colid+1)= 0;
            rowdata.at(colid+2)=-1;
            lpsolve::add_constraint(_component_lp.get(), &rowdata[0], ROWTYPE_LE,-target);
        }
    }

    void ComponentAssigner::build_integer_linear_program(){

        // lpsolve::delete_lp(_component_lp);
        // lpsolve::delete_lp(_component_lpcopy);
        _component_lp = unique_lp(lpsolve::make_lp(0,0), lpsolve::delete_lp);
        lpsolve::set_scaling(_component_lp.get(), SCALE_NONE);
        /* NOTE: *ALWAYS* use at least PRESOLVE_COLS. It allows for simple determination of
           possibly resonant systems, and significantly reduces the model size. */
        int presolvetype=PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP | PRESOLVE_MERGEROWS;
        presolvetype |= PRESOLVE_IMPLIEDSLK | PRESOLVE_REDUCEGCD | PRESOLVE_BOUNDS;
        //presolvetype=PRESOLVE_NONE;
        lpsolve::set_presolve(_component_lp.get(),presolvetype, lpsolve::get_presolveloops(_component_lp.get()) );

#if DEBUGPRINT==0
        lpsolve::set_verbose(_component_lp.get(),0);
#endif

        // Step 1) Add columns and set objective function to minimize
        set_atom_lonepair_penalties();
        set_bond_penalties();
        set_atom_charge_penalties();
        set_aromatic_ring_penalties();
        set_component_charge_penalty();

        /* The objective values for presolved columns get removed,
           so we copy the total objective vector and save it here */
        _component_objf.assign(lpsolve::get_Norig_columns(_component_lp.get())+1,0);
        lpsolve::get_row(_component_lp.get(),0,&_component_objf[0]);
        _component_objf[0]=0.0;

        // Step 2) Add rows (constraints/equalities)
        assert(lpsolve::set_add_rowmode(_component_lp.get(), true));

        add_indicator_constraints();
        add_atom_octet_and_charge_constraints();
        add_aromatic_ring_constraints();
        add_component_electron_constraint();

        assert(lpsolve::set_add_rowmode(_component_lp.get(), false));

        /* after calling lpsolve::solve when presolve is in effect, its impossible to change
           any of the constraints (eg total charge) and re-solve. Therefore, we copy the model
           and update the constraints before (re)solving the model */
        _component_lpcopy = unique_lp(lpsolve::copy_lp(_component_lp.get()), lpsolve::delete_lp);
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
        for (ComponentAssignerPtr ca : _component_assigners){
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
        for (Id cid : active){
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
        for (value_type const& cadata : solutions){
            cons.push_back(std::vector<double>());
            std::vector<double> &rowdata=cons.back();
            rowdata.resize(ncols+1,0);
            Id caidx=cadata.first;
            /* add new columns */
            for (entry_type const& soldata : cadata.second){
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
        for (std::vector<double> &rowdata : cons){
            lpsolve::add_constraint(qtotlp, &rowdata[0], ROWTYPE_EQ, 1);
        }
        lpsolve::set_add_rowmode(qtotlp, false);
        lpsolve::set_timeout(qtotlp, seconds_until_deadline());
        int status=lpsolve::solve(qtotlp);
        if ((status == SUBOPTIMAL) && (std::chrono::system_clock::now() > _deadline)) {
            status = TIMEOUT;
        }
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
            for (value_type const& cadata : solutions){
                Id caidx=cadata.first;
                int nset=0;
                for (entry_type const& soldata : cadata.second){
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
        if (status == TIMEOUT) {
            throw std::runtime_error("Unable to solve ILP. Timeout elapsed");
        }

        return _valid;
    }

    double BondOrderAssigner::getSolvedObjective(){
        if(!_valid){
            throw std::runtime_error("Cannot getSolvedObjective from invalid integer linear program solution."
                                     " Did you call solveIntegerLinearProgram first?");
        }
        double obj=0.0;
        for (ComponentAssignerPtr ca : _component_assigners){
            obj+= ca->getSolvedComponentObjective();
        }
        return obj;
    }

    void BondOrderAssigner::assignSolutionToAtoms(std::unordered_map<Id, std::unordered_map<Id, std::vector<int> > > &fc_groups,
                                                  std::unordered_map<Id, std::unordered_map<Id, std::vector<int> > > &bo_groups){

        if(!_valid){
            throw std::runtime_error("Cannot assignSolutionToAtoms with invalid integer linear program solution."
                                     " Did you call solveIntegerLinearProgram first?");
        }

        SolutionMap atominfo(_atominfo);
        SolutionMap bondinfo(_bondinfo);
        SolutionMap chargeinfo(_chargeinfo);

        /* Update presolved atominfo/bondinfo/chargeinfo with solved values from componentAssigners */
        for (ComponentAssignerPtr ca : _component_assigners){
            ca->extractComponentSolution(atominfo, bondinfo, chargeinfo);
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
        for (Id aid1 : _fragatoms){
            _mol->atom(aid1).formal_charge=0;
        }

        /* Assign bond orders and bond electron charge part here "- 0.5*Sum_j ( BondElectrons[ij] )" */
        auto qprop = BadId, oprop = BadId;
        if (compute_resonant_charge()) {
            qprop = _mol->addAtomProp("resonant_charge", FloatType);
            oprop = _mol->addBondProp("resonant_order", FloatType);
        }
        std::vector<double> resonant_charge(_mol->maxAtomId());

        std::unordered_map<Id, solutionResult> data =  bondinfo.extract();

        for (auto const& entry: data){
            Id bid = entry.first;
            bond_t & bond=_mol->bond(bid);
            Id aid1=bond.i;
            atom_t& atm1=_mol->atomFAST(aid1);
            Id aid2=bond.j;
            atom_t& atm2=_mol->atomFAST(aid2);

            /* Bond Orders */
            int order = entry.second.nonresonant;
            double resorder = entry.second.resonant;
            bond.order=order;
            if (compute_resonant_charge()) {
                _mol->bondPropValue(bid, oprop) = resorder;
            }

            /* Charges */
            atm1.formal_charge-=order;
            atm2.formal_charge-=order;
            resonant_charge[bond.i] -= resorder;
            resonant_charge[bond.j] -= resorder;
        }

        /* Take care of the "ValenceElectrons[i] - freeElectrons[i]" part of charge here */
        data = atominfo.extract();
        std::unordered_map<Id, solutionResult> qdata = chargeinfo.extract();
        for (auto const& entry: data){
            Id aid = entry.first;
            atom_t& atm=_mol->atom(aid);
            int nValence=DataForElement(atm.atomic_number).nValence;

            atm.formal_charge+= nValence - 2*entry.second.nonresonant;

            double resq = resonant_charge[aid];
            resq += nValence - 2*entry.second.resonant;
            if (compute_resonant_charge()) {
                _mol->atomPropValue(aid, qprop) = fabs(resq) < 1e-5 ? 0 : resq;
            }

            int iqtarget = 0;
            double dqtarget = 0.0;
            auto qresult = qdata.find(aid);
            if (qresult != qdata.end()){
                iqtarget=qresult->second.nonresonant;
                dqtarget=qresult->second.resonant;
            }
            if((atm.formal_charge != iqtarget) || fabs(resq - dqtarget)>1E-8 ){
                throw std::runtime_error("SolutionToAtoms failed. atm.formal_charge != iqtarget or fabs(resq - dqtarget)>1E-8");
            }
        }

        fc_groups = chargeinfo.generate_and_remove_resonance_groups();
        bo_groups = bondinfo.generate_and_remove_resonance_groups();

        /* Clean up {fc_groups:
            if the entry doesnt show up in bo_groups, make sure there is only one resonance and then
            move it to the default group
        */
        std::vector<Id> removed;
        for (auto const& kv: fc_groups){
            if ( kv.first == Id(-1) || bo_groups.count(kv.first) != 0) continue;
            /* This might create a default map if there wasnt one alread */
            auto & default_map = fc_groups[Id(-1)];
            for (auto const& kv2: kv.second){
                /* Pedantic: atom shouldnt be in the default map and should have a single resonance */
                assert(default_map.count(kv2.first) ==0 && kv2.second.size()==1);
               /* move the single atom, single resonance to the default list */
                default_map[kv2.first] = std::move(kv2.second);
            }
            removed.push_back(kv.first);
        }
        for (Id rgrp: removed){
            fc_groups.erase(rgrp);
        }



    }
}}
