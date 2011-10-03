#ifndef desres_msys_atomsel_keyword_predicate_hxx
#define desres_msys_atomsel_keyword_predicate_hxx

#include "keyword.hxx"
#include "predicate.hxx"

namespace desres { namespace msys { namespace atomsel {

  /* a predicate that evaluates true when corresponding values are
   * matched by at least one of its literals or ranges. */
  class KeywordPredicate : public Predicate {
    KeywordPtr key;
    std::set<Literal>   literals;
    std::set<Range>     ranges;
    std::vector<Regex>  regexes;

    public:
    KeywordPredicate(KeywordPtr _key) : key(_key) {}

    void eval( Selection& s ) { key->select(s,literals,ranges,regexes); }
    void dump(std::ostream& str) const;

    void addLiteral( const Literal& lit ) { literals.insert(lit); }
    void addRange( const Range& range) { ranges.insert(range); }
    void addRegex( const Regex& regex ) { regexes.push_back(regex); }
  };

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
