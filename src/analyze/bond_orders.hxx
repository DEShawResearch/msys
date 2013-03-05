#ifndef desres_msys_analyze_bond_orders_hxx
#define desres_msys_analyze_bond_orders_hxx

#include "../system.hxx"
#include "bondFilters.hxx"
#include <boost/noncopyable.hpp>


/* forward declare so we dont pollute the rest of our code with lpsolves excessive #defines */
namespace lpsolve{
    class _lprec;
}

namespace desres { namespace msys {

    struct electronRange {
        int lb;
        int ub;
    };
    typedef std::map<Id, electronRange> electronMap;

    struct solutionValues {
        solutionValues():nonresonant(0),resonant(0){};
        solutionValues(int i,double d):nonresonant(i),resonant(d){};
        int nonresonant;
        double resonant;
    };
    typedef std::map<Id, solutionValues> solutionMap;

    class ComponentAssigner;
    typedef boost::shared_ptr<ComponentAssigner> ComponentAssignerPtr;


    class BondOrderAssigner : public boost::noncopyable {
        friend class ComponentAssigner;

#define GETSETREBUILDVAR(vtype, name)                                   \
    private:                                                            \
        vtype name;                                                     \
    public:                                                             \
        const vtype& get_##name() const { return name; }               \
        void set_##name(const vtype& newval) { _needRebuild=true; _valid=false; name = newval; }

        /* Weighting factor for the 2 atomic electronegativities in a bond. 
           Set to 0.5 for equal weights, 0.0 for no contribution from lower eneg atom*/
        GETSETREBUILDVAR(double, bond_eneg_factor);
        /* scale factor for charged atom penalty. 
           Set to 0.0 to remove penalty for forming charged atoms  */
        GETSETREBUILDVAR(double, atom_charge_penalty_factor);
        /* scale factor for charged component penalty.  
           Set to 0.0 to remove penalty for forming charged components  */
        GETSETREBUILDVAR(double, component_charge_penalty_factor);
        /* scale factor for not forming aromatic rings. 
           Set to 0.0 to remove penalty for not forming aromatic rings */
        GETSETREBUILDVAR(double, ring_penalty_factor);
        /* largest allowed absolute atom charge (if gen_charge_penalty_for_atom) */
        GETSETREBUILDVAR(unsigned, max_atom_charge);
        /* largest allowed absolute component (resonant system) charge */
        GETSETREBUILDVAR(unsigned, max_component_charge);
        /* Rings with planarity values < planarity_tolerance are considered planar */
        GETSETREBUILDVAR(double, planarity_tolerance);
        /* Penalty for forming a hextet atom. Should be large so hextets are formed as a last resort */
        GETSETREBUILDVAR(double, atom_hextet_penalty);
#undef GETSETREBUILDVAR
    public:
        BondOrderAssigner(){
            _needRebuild=true; 
            _valid=false;
            planarity_tolerance=0.1;
            /* golden ratio less contribution from lower eneg atom */
            bond_eneg_factor=0.3819660112501059;
            atom_charge_penalty_factor=1.0;
            component_charge_penalty_factor=1.0;
            ring_penalty_factor=1.0;
            atom_hextet_penalty=50.0;
            max_atom_charge=2;
            max_component_charge=6; 
            _filter=NULL;
        };
        ~BondOrderAssigner(){ delete _filter;};


        static boost::shared_ptr<BondOrderAssigner> create(SystemPtr sys, 
                                                           IdList const& fragments);

        bool allow_hextet_for_atom(Id aid1);
        void reset();

        void setTotalCharge(int qTotal);
        void unsetTotalCharge();
        int getSolvedTotalCharge();

        void setComponentCharge(Id component, int qTotal);
        void unsetComponentCharge(Id component);
        int  getSolvedComponentCharge(Id component);

        bool solveIntegerLinearProgram();
        void assignSolutionToAtoms();

        int getTotalValence(){return _totalValence;}
        double getSolvedObjective();

    protected:
        int max_free_pairs(const Id aid);
        int max_bond_order(const Id aid);
      
        void presolve_octets(IdList& unsolved);
        void rebuild();

        /* setup during create() */
        bool _needRebuild;      // rebuild the rings and component assigners before solving ?
        bool _valid;            // Do we currently have a valid solution? 
        int _totalValence;      // Total valence elctron count
        int _total_charge;      // total charge of fragment
        int _presolved_charge;  // charge of presolved components
        bool _total_charge_set; // did we set the total charge?
        SystemPtr _mol; 
        bondFilter *_filter;

        IdList _fragatoms;
        /* These contain the presolved atom lp, bond orders and atomic charges */
        solutionMap _atominfo;
        solutionMap _bondinfo;
        solutionMap _chargeinfo;
        electronMap _atom_lp;
        electronMap _bond_order;

        /* filled during rebuild */
        typedef std::pair<double, IdList> RingPair;
        typedef std::vector<RingPair> RingList;
        RingList _planar_rings;
        RingList _nonplanar_rings;
        std::vector<ComponentAssignerPtr> _component_assigners;
        std::map<Id, int> _fixed_component_charges;

    };
    typedef boost::shared_ptr<BondOrderAssigner> BondOrderAssignerPtr;
    typedef boost::weak_ptr<BondOrderAssigner> BondOrderAssignerWeakPtr;


    class ComponentAssigner : public boost::noncopyable {
        static const int MIN_INVALID_ILP_COL=0;
        ComponentAssigner(): _component_lp(NULL), _component_lpcopy(NULL){};
    public:

        ~ComponentAssigner();
        
        static boost::shared_ptr<ComponentAssigner> create(BondOrderAssignerPtr boa,
                                                           IdList const& comp, 
                                                           Id cid );
        /* lpsolve helper routines */
        static int add_column_to_ilp(lpsolve::_lprec *lp, std::string const& colname, 
                                     double penalty, double lb, double ub);
        static void get_ilp_solution(lpsolve::_lprec *lp, std::vector<int> &solution);

        void build_integer_linear_program();

        void reset();

        void setComponentCharge(int qTotal);
        void unsetComponentCharge();
        int  getSolvedComponentCharge();

        bool solveComponentIntegerLinearProgram();
        void extractComponentSolution(solutionMap &atominfo,
                             solutionMap &bondinfo, 
                             solutionMap &chargeinfo);

        int getComponentValence(){return _component_valence_count;}
        double getSolvedComponentObjective();

    protected:
        struct ringAtomInfo{
            ringAtomInfo(): nBonds(0), lonePairIdx(0),
                            bondIdxPrevious(0), bondIdxNext(0),
                            bondIdxExoCyclic(0){};
            Id aid;
            int nBonds;
            int lonePairIdx;
            int bondIdxPrevious;
            int bondIdxNext;
            int bondIdxExoCyclic;
        };
        struct resPathInfo{
            resPathInfo(): atom0_daidx(0), atom0_lpcol(0), atom0_qcol(0),
                           atom1_daidx(0), atom1_lpcol(0), atom1_qcol(0){};
            int atom0_daidx; 
            int atom0_lpcol;
            int atom0_qcol;
            int atom1_daidx;
            int atom1_lpcol;
            int atom1_qcol;
            std::vector< std::pair<int,int> > bondColsMaxOrders;
        };
        

        /* This has to be a weak_ptr otherwise we can never 
         * destroy parent if we have any componentAssigners */
        BondOrderAssignerWeakPtr _parent;
        Id _component_id;
        std::set<Id> _component_atoms_present;
        
        int _component_valence_count;
        // bool _component_has_expanded_octets;
        
        typedef std::map<Id, int> ilpMap;
        ilpMap _component_atom_cols;
        ilpMap _component_atom_charge_cols;
        ilpMap _component_atom_hextet_cols;

        ilpMap _component_bond_cols;
        ilpMap _component_ring_cols;
        int _component_charge_col;
        int _component_charge_row;

        bool _component_solution_valid;
        std::vector<double> _component_objf;
        std::vector<int> _component_solution;
        std::vector<double> _component_resonant_solution;
        
        lpsolve::_lprec *_component_lp;
        lpsolve::_lprec *_component_lpcopy;
        
        void set_atom_hextet_penalties();

        bool gen_charge_penalty_for_atom(Id aid1);
        void set_atom_charge_penalties();

        void set_atom_lonepair_penalties();
        void set_bond_penalties();
        void set_aromatic_ring_penalties();
        void set_component_charge_penalty();
        
        void break_ring_symmetry();
        
        void add_atom_octet_and_charge_constraints();
        void generate_ring_constraint( Id ridx, int rcolid, double &target, 
                                       std::vector<double> &rowdata);
        void add_aromatic_ring_constraints();
        void add_component_electron_constraint();
        
        /* Resonance Generation */
        void check_resonance_path(resPathInfo const& respath,
                                  int direction, bool xferLonePair,
                                  std::vector<int> const& soln,
                                  std::vector<int> & newsoln);
        
        unsigned check_aromatic_path(std::vector<ringAtomInfo> const& ringData,
                                std::vector<int> const& soln,
                                std::vector<int> & newsoln);
        
        void initialize_donor_acceptors(std::vector<unsigned> & da_state,  
                                       std::vector<int> & transfertype,
                                       std::vector<resPathInfo> & resonance_paths );
        
        void initialize_aromatic_ringdata(std::vector< std::vector< ringAtomInfo > > &rings_to_process);
        
        
        void update_aring_constraint_for_resonance_form(std::vector<std::vector<double> > const& arom_cons,
                                                        std::vector<int> const& old_soln, 
                                                        std::vector<int> & new_soln);
        
        void generate_resonance_forms_to_check(std::vector<std::vector<double> > const& arom_cons,
                                               std::set< std::vector<int> > &rforms);
        
        void sum_resonance_solutions(std::set< std::vector<int> > const& rforms);
        
        void extract(ilpMap const& im, bool isCharge, solutionMap & sm);
        
    };

}}
#define DEBUGPRINT  0
#define DEBUGPRINT1 0
#define DEBUGPRINT2 0
#endif
