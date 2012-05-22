#ifndef desres_msys_override_hxx
#define desres_msys_override_hxx

#include "param_table.hxx"

namespace desres { namespace msys {

    class OverrideTable;
    typedef boost::shared_ptr<OverrideTable> OverrideTablePtr;

    class OverrideTable {

        /* which parameter table's entries are being overridden */
        ParamTablePtr   _target;

        /* storage for override values */
        ParamTablePtr   _params;

        /* all the overrides.  Keys have first <= second */
        typedef std::map<IdPair, Id> OverrideMap;
        OverrideMap     _map;

        /* constructor: the parameter table we override */
        explicit OverrideTable( ParamTablePtr target );

    public:
        /* create an override table */
        static OverrideTablePtr create(ParamTablePtr target);

        /* clear map and remove all references to target */
        void clear();

        /* target */
        ParamTablePtr target() const { return _target; }
        
        /* override parameters */
        ParamTablePtr params() const { return _params; }

        /* get parameter for the given pair of ids in target;
         * BadId if not present.  The input arguments will be sorted,
         * so the order in which they are provided is irrelevant */
        Id get(IdPair params) const;

        /* remove the given override if present */
        void del(IdPair params);

        /* add or replace an override */
        void set(IdPair params, Id param);

        /* number of distinct overrides */
        Id count() const;

        /* list of all the overrides */
        std::vector<IdPair> list() const;
    };

}}

#endif
