#ifndef desres_msys_term_table_hxx
#define desres_msys_term_table_hxx

#include "override.hxx"

namespace desres { namespace msys {

    /* The TermTable class holds forcefield terms, which are a mapping 
     * from n-tuples of atoms to a parameter set.  */
    
    class System;
    typedef std::shared_ptr<System> SystemPtr;

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

    class TermTable : public std::enable_shared_from_this<TermTable> {
        std::weak_ptr<System> _system;
        ParamTablePtr   _params;
    
        /* number of atoms in each term */
        Id      _natoms;

        /* number of dead terms */
        Id      _ndead;
    
        /* all the terms. */
        IdList  _terms;

        /* terms are live if their first atom is not BadId */
        bool _alive(Id term) const {
            return BadId != _terms.at(term*(1+_natoms));
        }

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

        std::map<String, String> _unused_keep_for_binary_compatibility;

    public:
        template <class Archive>
        void serialize(Archive & archive) {
            archive(_system, _params, _natoms, _ndead, _terms, _props, _overrides, category);
        }

        // default constructor for serialization only
        TermTable() {}

        TermTable( SystemPtr system, Id natoms, 
                   ParamTablePtr ptr = ParamTablePtr() );

        /* delete all terms and remove from parent system */
        void destroy();

        SystemPtr system() const { return _system.lock(); }
        ParamTablePtr params() { return _params; }
        inline ParamTablePtr props() { return _props; }
        Id atomCount() const { return _natoms; }

        /* category describing what sort of TermTable we have */
        Category category;

        /* name of this table in the parent system */
        String name() const;

        /* rename table */
        void rename(String const& name);

        /* reserve space */
        void reserve(size_t nterms) { _terms.reserve(nterms * (_natoms+1)); }

        class term_t {
        protected:
            const Id*   _ptr;
            Id          _size;
            Id          _id;

            term_t(const Id* ptr, Id size, Id id) 
            : _ptr(ptr), _size(size), _id(id) {}

        public:
            term_t() 
            : _ptr(), _size(), _id(BadId) {}

            const Id* atoms() const { return _ptr; }
            Id const& atom(int i) const { return _ptr[i]; }
            Id size() const { return _size; }
            Id id() const { return _id; }
            Id param() const { return _ptr[_size]; }
        };

        /* const iterator access over all terms (even the dead ones) */
        class const_iterator : protected term_t {
            friend class TermTable;

            const_iterator(const Id* ptr, Id size, Id id)
            : term_t(ptr+id*(size+1), size, id) {}

            void increment() { 
                _ptr += _size+1;
                ++_id;
            }

            bool equal(const_iterator const& c) const { return _ptr==c._ptr; }
            const term_t& dereference() const { return *this; }

        public:
            const_iterator() : term_t() {}
            const term_t& operator*() const { return dereference(); }
            const term_t* operator->() const { return &dereference(); }
            const_iterator& operator++() { increment(); return *this; }
            bool operator==(const_iterator const& c) const { return equal(c); }
            bool operator!=(const_iterator const& c) const { return !equal(c); }
        };

        /* iterator pointing to first element */
        const_iterator begin() const {
            if (_terms.empty()) return const_iterator();
            return const_iterator(&_terms[0], _natoms, 0);
        }
        /* iterator pointing past last element */
        const_iterator end() const {
            if (_terms.empty()) return const_iterator();
            return const_iterator(&_terms[0], _natoms, maxTermId());
        }

        /* Operations on the set of terms */
        IdList terms() const;
        Id termCount() const {
            return maxTermId() - _ndead;
        }
        Id maxTermId() const {
            return _terms.size()/(1+_natoms);
        }

        bool hasTerm(Id term) const {
            return term<maxTermId() && _alive(term);
        }

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

        /* Return the ids of terms which contain only the given atoms;
         * these are the terms which would be included in a clone().
         *
         * Internally, this just calls findWithAny(), then fiters the result.
         */
        IdList findWithOnly(IdList const& ids);

        /* return the ids of the terms containing _exactly_ the given ids,
         * in the given order.
         *
         * Internally, this just calls findWithAll(), then filters the result.
         */
        IdList findExact(IdList const& ids);

        /* Operations on individual terms */
        Id param(Id term) const;
        void setParam(Id term, Id param);
        IdList atoms(Id term) const;
        Id atom(Id term, Id index) const;
        const Id* atomsFAST(Id term) const {
            return &_terms[term*(1+_natoms)];
        }
        Id paramFAST(Id term) const {
            return _terms[(1+term)*(1+_natoms)-1];
        }

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
        OverrideTablePtr overrides() { return _overrides; }

        /* replace the parameter table with a new one.  The param of all
         * terms will be set to BadId, even if the new ParamTable is the
         * same as the old one.  The override table will also be
         * replaced with a new one pointing to the new ParamTable. */
        void resetParams(ParamTablePtr params);

        inline ParamTablePtr termProps() { return _props; }
    };

    typedef std::shared_ptr<TermTable> TermTablePtr;

    // Construct a new table with terms sorted by atom ids.
    // Replace the table in the parent system with the sorted version.
    // Return the new table.
    TermTablePtr ReplaceTableWithSortedTerms(TermTablePtr table);
}}

#endif

