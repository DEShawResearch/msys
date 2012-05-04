#ifndef desres_msys_tuples_hxx
#define desres_msys_tuples_hxx

#include "term_table.hxx"

namespace desres { namespace msys {

    /* Inputs: base -- a fan-1 term table 
     *         combined - a table with columns param1, param2, ... paramN
     *                    as well as a subset of the columns in base.params.
     *          tuples -- a fan-N term table sharing a param table with base.
     *
     * Side effects: 
     *         For each param1, param2, ... paramN tuple in combined,
     *         add to tuples all possible atom tuples atom1, atom2, ... atomN 
     *         where atomK is drawn from the set of all atoms with
     *         param paramK.  New parameter entries will be added to base.
     */
    void CreateTuplesFromCombined( TermTablePtr base,
                                   ParamTablePtr combined,
                                   TermTablePtr tuples );

    /* Attempt to invert the pairs creation.  This can be done iff
     * all tuples in pairs with the same param assignment in base have
     * the equivalent param assignment in pairs. 
     *
     * Inputs: base - a fan-1 term table.
     *         tuples - a fan-N term table sharing a param table with base.
     *
     * Output: a ParamTable with the same form as above.
     *
     * Side effects: none.
     */
    ParamTablePtr CreateCombinedFromTuples( TermTablePtr base,
                                            TermTablePtr tuples );
}}

#endif
