#ifndef desres_msys_atomsel_predicate_hxx
#define desres_msys_atomsel_predicate_hxx

#include "selection.hxx"
#include <ostream>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace desres { namespace msys { namespace atomsel {

  /* a predicate ANDs the items in the input selection with a condition
   * in a particular context. */
  class Predicate : public boost::enable_shared_from_this<Predicate> {
  public:
    virtual ~Predicate() {}
    virtual void eval( Selection& s ) = 0;
    virtual void dump(std::ostream& str) const = 0;
  };

  typedef boost::shared_ptr<Predicate> PredicatePtr;

  class StringPredicate : public Predicate {
    public:
        virtual void add( std::string const& elem ) = 0;
  };

  /* All flags true */
  PredicatePtr all_predicate();

  /* All flags false */
  PredicatePtr none_predicate();

  /* the AND predicate: S -> S AND (L AND R ) */
  PredicatePtr and_predicate( PredicatePtr L, PredicatePtr R );

  /* The OR predicate: S -> S AND (L OR R) */
  PredicatePtr or_predicate( PredicatePtr L, PredicatePtr R );

  /* the NOT predicate: S -> S AND ( NOT X ) */
  PredicatePtr not_predicate( PredicatePtr X );

}}} // ns

#endif
