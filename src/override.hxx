#ifndef desres_msys_override_hxx
#define desres_msys_override_hxx

#include "term_table.hxx"

namespace desres { namespace msys {

    /* Inputs: base   -- a single atom table with one term for every atom.
     *         tuples -- a two-atom table with distinct entries
     *
     * Output: a mapping from the overridden base parameter ids to the
     *         id of their overridden value in tuples.
     *         It is an error if terms in tuples with the same param ids
     *         in base have different parameter values in tuples.
     *
     * Side effects: none.
     */

    typedef std::map<std::pair<Id,Id>, Id> OverrideMap;
    OverrideMap FindOverridesFromTuples( TermTablePtr base,
                                         TermTablePtr tuples );

    /* Add terms to override table tuples based on given overrides.
     * The key ids in overrides must point to existing params in base;
     * the value id points to an existing param in tuples. */
    void MakeTuplesFromOverrides( OverrideMap const& overrides,
                                  TermTablePtr base,
                                  TermTablePtr tuples );
}}

#endif
