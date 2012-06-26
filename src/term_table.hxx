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

    class OverrideTable;
    typedef boost::shared_ptr<OverrideTable> OverrideTablePtr;

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

        /* overrides for the parameters of this table */
        OverrideTablePtr _overrides;

        /* mapping from atom id to list of terms having that id. */
        std::vector<IdList> _index;
        /* max term id from last update */
        Id _maxIndexId;

        /* the index is updated only by the find() operations, and by 
         * delTerm if an index has already been created.  */
        void update_index();

    public:
        TermTable( SystemPtr system, Id natoms, 
                   ParamTablePtr ptr = ParamTablePtr() );

        /* delete all terms and remove from parent system */
        void destroy();

        SystemPtr system() const { return _system.lock(); }
        ParamTablePtr params() { return _params; }
        Id atomCount() const { return _natoms; }

        /* category describing what sort of TermTable we have */
        Category category;

        /* name of this table in the parent system */
        String name() const;

        /* rename table */
        void rename(String const& name);
    
        /* Operations on the set of terms */
        IdList terms() const;
        Id termCount() const;
        Id maxTermId() const;
        bool hasTerm(Id term) const;
        Id addTerm(const IdList& atoms, Id param);
        void delTerm(Id id);

        /* delete all terms t containing atom id atm i the atoms list.  */
        void delTermsWithAtom(Id atm);

        /* return the ids of the terms which contain _all_ of the
         * given atoms, in any order.  */
        IdList findWithAll(IdList const& ids);

        /* Return the ids of the terms which contain _any_ (at least one)
         * of the given atoms, in any order. */
        IdList findWithAny(IdList const& ids);

        /* return the ids of the terms containing _exactly_ the given ids,
         * in the given order.
         *
         * Internally, this just calls findAllIds(), then filters the result.
         */
        IdList findExact(IdList const& ids);

        /* Operations on individual terms */
        Id param(Id term) const;
        void setParam(Id term, Id param);
        IdList atoms(Id term) const;
        Id atom(Id term, Id index) const;

        /* look up the value of a property of the term from the associated
         * ParamTable */
        ValueRef propValue(Id term, Id index);
        ValueRef propValue(Id term, String const& name);

        /* Operations on term properties */
        Id termPropCount() const;
        String termPropName(Id i) const;
        ValueType termPropType(Id i) const;
        Id termPropIndex(String const& name) const;
        Id addTermProp(String const& name, ValueType type);
        void delTermProp(Id index);
        ValueRef termPropValue(Id term, Id index);
        ValueRef termPropValue(Id term, String const& name);

        /* reassign param to a member of the set of distinct parameters. */
        void coalesce();

        /* get the override table */
        OverrideTablePtr overrides() const { return _overrides; }
    };

    typedef boost::shared_ptr<TermTable> TermTablePtr;

}}

#endif

