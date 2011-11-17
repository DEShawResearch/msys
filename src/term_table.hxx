#ifndef desres_msys_term_table_hxx
#define desres_msys_term_table_hxx

#include "types.hxx"
#include "param_table.hxx"
#include <boost/enable_shared_from_this.hpp>

namespace desres { namespace msys {

    /* The TermTable class holds forcefield terms, which are a mapping 
     * from n-tuples of atoms to a parameter set.  */
    
    class System;
    typedef boost::shared_ptr<System> SystemPtr;

    /* Every term table should be assigned a category based on its type. */
    enum Category {
        NO_CATEGORY = 0,
        BOND        = 1,
        CONSTRAINT  = 2,
        VIRTUAL     = 3,
        POLAR       = 4,
        NONBONDED   = 5,
        EXCLUSION   = 6
    };

    /* convert category to string */
    std::string print(Category const& c);

    /* convert string to category */
    Category parse(std::string const& s);

    class TermTable : public boost::enable_shared_from_this<TermTable> {
        boost::weak_ptr<System> _system;
        ParamTablePtr   _params;
    
        /* number of atoms in each term */
        Id          _natoms;
    
        /* all the terms. */
        typedef std::vector<Id> TermList;
        TermList    _terms;
        IdSet   _deadterms; 

        /* extra term properties, analogous to extra atom properties.
         * These are currently used for two things in dms files: 
         * the "constrained" field to mark constrained stretch and 
         * angle terms, and for the position of position restraints. */
        ParamTablePtr  _props;

        /* parameter id for alchemical state.  We expect most terms,
         * and indeed most term tables, will not be alchemical, so 
         * we store these ids in a map.  We never put BadId into the
         * map; thus, we can return BadId for a term's alchemical
         * parameter to indicate no alchemical state exists for that
         * term, and the size of the map can be used to indicate whether
         * any part of the table is alchemical. */
        typedef std::map<Id,Id> ParamMap;
        ParamMap    _paramB;
    
    public:
        TermTable( SystemPtr system, Id natoms, 
                   ParamTablePtr ptr = ParamTablePtr() );
        SystemPtr system() const { return _system.lock(); }
        ParamTablePtr params() { return _params; }
        Id atomCount() const { return _natoms; }

        /* true if there are any alchemical terms */
        bool alchemical() const;

        /* category describing what sort of TermTable we have */
        Category category;

        /* name of this table in the parent system */
        String name() const;
    
        /* Operations on the set of terms */
        IdList terms() const;
        Id termCount() const;
        Id maxTermId() const;
        bool hasTerm(Id term) const;
        Id addTerm(const IdList& atoms, Id param);
        void delTerm(Id id);

        /* delete all terms t containing atom id atm i the atoms list. 
         * This operation requires a full scan of the entire set of
         * terms, so use with care.  */
        void delTermsWithAtom(Id atm);
    
        /* Operations on individual terms */
        Id param(Id term) const;
        Id paramB(Id term) const;
        void setParam(Id term, Id param);
        void setParamB(Id term, Id param);
        IdList atoms(Id term) const;
        Id atom(Id term, Id index) const;

        /* Operations on term properties */
        Id termPropCount() const;
        String termPropName(Id i) const;
        ValueType termPropType(Id i) const;
        Id termPropIndex(String const& name) const;
        Id addTermProp(String const& name, ValueType type);
        void delTermProp(Id index);
        ValueRef termPropValue(Id term, Id index);
        ValueRef termPropValue(Id term, String const& name);
    };

    typedef boost::shared_ptr<TermTable> TermTablePtr;

}}

#endif

