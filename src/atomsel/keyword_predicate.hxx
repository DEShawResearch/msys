#ifndef desres_msys_atomsel_keyword_predicate_hxx
#define desres_msys_atomsel_keyword_predicate_hxx

#include "keyword.hxx"
#include "predicate.hxx"

namespace desres { namespace msys { namespace atomsel {

    /* Select from list of targets */
    PredicatePtr keyword_predicate(KeywordPtr key, TargetList* targets);

    /* a boolean predicate is just a keyword predicate with the single 
     * literal value of "1". */
    PredicatePtr boolean_predicate( KeywordPtr key );

    /* A predicate which is true for elements whose keyword value is in the
     * union of values from a subselection */
    PredicatePtr same_keyword_predicate( KeywordPtr key, 
                                         const Selection& all,
                                         PredicatePtr sub );

}}}

#endif
