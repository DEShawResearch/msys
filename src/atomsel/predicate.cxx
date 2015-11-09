#include "predicate.hxx"

namespace {
  using namespace desres::msys::atomsel;

  struct AllPredicate : Predicate {
    void eval(Selection& s) {}
  };

  struct NonePredicate : Predicate {
    void eval(Selection& s) { s.clear(); }
  };

  struct AndPredicate : Predicate {
    PredicatePtr left, right;
    AndPredicate( PredicatePtr A, PredicatePtr B ) : left(A), right(B) {}
    void eval(Selection& s) {
      left->eval(s);
      right->eval(s);
    }
  };
  struct OrPredicate : Predicate {
    PredicatePtr left, right;
    OrPredicate( PredicatePtr A, PredicatePtr B ) : left(A), right(B) {}
    void eval(Selection& s) {
      Selection s2(s);
      left->eval(s);
      right->eval(s2);
      s.add(s2);
    }
  };
  struct NotPredicate : Predicate {
    PredicatePtr sub;
    NotPredicate( PredicatePtr A ) : sub(A) {}
    void eval(Selection& s) {
      Selection s2(s);
      sub->eval(s2);
      s.subtract(s2);
    }
  };
}

namespace desres { namespace msys { namespace atomsel {

  PredicatePtr all_predicate() {
      static PredicatePtr p(new AllPredicate);
      return p;
  }

  PredicatePtr none_predicate() {
      static PredicatePtr p(new NonePredicate);
      return p;
  }

  PredicatePtr and_predicate( PredicatePtr L, PredicatePtr R ) {
      if (!L || !R) return PredicatePtr();
      return PredicatePtr(new AndPredicate(L,R));
  }

  PredicatePtr or_predicate( PredicatePtr L, PredicatePtr R ) {
      if (!L || !R) return PredicatePtr();
      return PredicatePtr(new OrPredicate(L,R));
  }

  PredicatePtr not_predicate( PredicatePtr X ) {
      if (!X) return PredicatePtr();
      return PredicatePtr(new NotPredicate(X));
  }

}}}
