#ifndef desres_msys_analyze_bond_orders_hxx
#define desres_msys_analyze_bond_orders_hxx

#include "../system.hxx"
#include "bondFilters.hxx"
#include <boost/noncopyable.hpp>
#include <array>

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
        /* scale for atoms lone pair electronegativity. 
           Set to 0.0 to remove penalty for forming positivly charged atoms  */
        GETSETREBUILDVAR(double, atom_lone_pair_scale);
        /* penalty for positivly charged atoms. 
           Set to 0.0 to remove penalty for forming positivly charged atoms  */
        GETSETREBUILDVAR(double, atom_plus_charge_penalty);
        /* penalty for negativly charged atoms. 
           Set to 0.0 to remove penalty for forming negativly charged atoms  */
        GETSETREBUILDVAR(double, atom_minus_charge_penalty);
        /* scale factor for charged atoms. 
           Set to 0.0 to remove scale factor forming  charged atoms  */
        GETSETREBUILDVAR(double, atom_charge_scale);
        /* penalty for charged terminal atoms. 
           Set to 0.0 to remove penalty for forming charged terminal atoms
           without adjacent compensating charge */
        GETSETREBUILDVAR(double, atom_terminal_charge_penalty);
        /* penalty for removing bond to "hypervalent" atoms. 
           Set to 0.0 to remove penalty for removing bonds to hypervalent atoms */
        GETSETREBUILDVAR(double, bond_fisure_penalty);
        /* penalty for positivly charged components.  
           Set to 0.0 to remove penalty for forming positivly charged components  */
        GETSETREBUILDVAR(double, component_plus_charge_penalty);
        /* penalty for negativly charged components.  
           Set to 0.0 to remove penalty for forming negativly charged components  */
        GETSETREBUILDVAR(double, component_minus_charge_penalty);
        /* penalty for not forming aromatic rings. 
           Set to 0.0 to remove penalty for not forming aromatic rings */
        GETSETREBUILDVAR(double, aromatic_ring_penalty);


        /* largest allowed absolute atom charge (if gen_charge_penalty_for_atom) */
        GETSETREBUILDVAR(unsigned, absmax_atom_charge);
        /* largest allowed absolute component (resonant system) charge */
        GETSETREBUILDVAR(unsigned, max_component_charge);

#undef GETSETREBUILDVAR
    public:
        BondOrderAssigner(){
            _needRebuild=true; 
            _valid=false;
            atom_lone_pair_scale=0.5;
            atom_plus_charge_penalty=0.25;
            atom_minus_charge_penalty=0.50;
            atom_charge_scale=1.50;
            atom_terminal_charge_penalty=1.0;
            bond_fisure_penalty=2.0;
            component_plus_charge_penalty=0.75;
            component_minus_charge_penalty=0.25;
            aromatic_ring_penalty=0.25;
            absmax_atom_charge=2;
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
        std::array<int,3> minmax_bond_order(const Id aid0, const Id aid1);
      
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

        MultiIdList _rings;
        std::set<Id> _ringAtoms;

        /* filled during rebuild */
        std::vector<ComponentAssignerPtr> _component_assigners;
        std::map<Id, int> _fixed_component_charges;

    };
    typedef boost::shared_ptr<BondOrderAssigner> BondOrderAssignerPtr;
    typedef boost::weak_ptr<BondOrderAssigner> BondOrderAssignerWeakPtr;


    class ComponentAssigner : public boost::noncopyable {
        static const int MIN_INVALID_ILP_COL=0;
        ComponentAssigner(): _component_lp(NULL), _component_lpcopy(NULL), _component_reslp(NULL){};
    public:

        ~ComponentAssigner();
        
        static boost::shared_ptr<ComponentAssigner> create(BondOrderAssignerPtr boa,
                                                           IdList const& comp, 
                                                           Id cid );
        /* lpsolve helper routines */
        static int add_column_to_ilp(lpsolve::_lprec *lp, std::string const& colname, 
                                     double penalty, double lb, double ub);
        static void get_ilp_solution(lpsolve::_lprec *lp, std::vector<int> &solution);
        double get_ilp_objective(std::vector<int> const& solution);

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
        struct ilpAtom {
            int ilpCol;
            int ilpLB;
            int ilpUB;
            int qCol;
            int qCtr;
            std::set<int> qTerm;
            ilpAtom(){};
            ilpAtom(int c,int l, int u):ilpCol(c),ilpLB(l),ilpUB(u),
                                        qCol(0),qCtr(0){};
        };
        typedef std::map<Id, ilpAtom> ilpAtomMap;
        
        struct ilpBond {
            int ilpCol;
            int ilpLB;
            int ilpUB;
            ilpBond(){};
            ilpBond(int c,int l, int u):ilpCol(c),ilpLB(l),ilpUB(u){};
        };
        typedef std::map<Id, ilpBond> ilpBondMap;
        typedef std::map<Id,int> ilpRingMap;

        /* This has to be a weak_ptr otherwise we can never 
         * destroy parent if we have any componentAssigners */
        BondOrderAssignerWeakPtr _parent;
        Id _component_id;
        std::set<Id> _component_atoms_present;
        
        int _component_valence_count;
        
        ilpAtomMap _component_atom_cols;
        ilpBondMap _component_bond_cols;
        ilpRingMap _component_ring_cols;

        int _component_charge_col;
        int _component_charge_row;

        bool _component_solution_valid;
        std::vector<double> _component_objf;
        std::vector<int> _component_solution;
        std::vector<double> _component_resonant_solution;
        
        lpsolve::_lprec *_component_lp;
        lpsolve::_lprec *_component_lpcopy;
        lpsolve::_lprec *_component_reslp;      

        void set_atom_charge_penalties();
        void set_atom_terminal_charge_penalties();
        void set_atom_lonepair_penalties();
        void set_bond_penalties();
        void set_aromatic_ring_penalties();
        void set_component_charge_penalty();
        
        void add_charge_vector_for_atom(BondOrderAssignerPtr parent,
                                        Id aid, std::vector<double> &row);        void add_indicator_constraints();
        void add_atom_octet_and_charge_constraints();
        void add_aromatic_ring_constraints();
        void add_component_electron_constraint();
        
        /* Resonance Generation */
        void generate_resonance_forms();
          
    };

}}
#define DEBUGPRINT  0
#define DEBUGPRINT1 0
#define DEBUGPRINT2 0
#endif
