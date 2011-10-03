#include "keyword_predicate.hxx"
#include <stdexcept>

namespace {
  using namespace desres::msys::atomsel;
  class SamePredicate : public Predicate {
    KeywordPtr key;
    Selection all;
    PredicatePtr sub;
    public:
    SamePredicate( KeywordPtr _key, const Selection& _all, PredicatePtr _sub)
      : key(_key), all(_all), sub(_sub) {}

    void eval(Selection& s);
    void dump(std::ostream& str) const {
      str << "same " << key->name << " as ";
      sub->dump(str);
    }
  };
}

template <typename T>
static void select( const std::vector<T>& vals,
    const Selection& x,
    Selection& s ) {

  std::set<T> u;
  int i,n=x.size();
  for (i=0; i<n; i++) if (x[i]) u.insert(vals[i]);
  for (i=0; i<n; i++) s[i] &= u.count(vals[i]);
}

void SamePredicate::eval( Selection& s ) {
  /* evaluate full subselection */
  Selection x(all);
  sub->eval(x);

  /* we'll need values from the union of x and s */
  Selection xs(x);
  xs.add(s);

  if (key->type==KEY_INT) {
    std::vector<Int> vals(xs.size());
    key->iget(xs,vals);
    select( vals, x, s );
  } else if (key->type==KEY_DBL) {
    std::vector<Dbl> vals(xs.size());
    key->dget(xs,vals);
    select( vals, x, s );
  } else if (key->type==KEY_STR) {
    std::vector<Str> vals(xs.size());
    key->sget(xs,vals);
    select( vals, x, s );
  } else {
    throw std::runtime_error("bad keyword type");
  }
}

namespace desres { namespace msys { namespace atomsel {

  void KeywordPredicate::dump(std::ostream& str) const {
    str << key->name << " ";
    for (std::set<Literal>::const_iterator i=literals.begin();
        i!=literals.end(); ++i) {
      str << *i << " ";
    }
    for (std::set<Range>::const_iterator i=ranges.begin();
        i!=ranges.end(); ++i) {
      str << i->first << " to " << i->second << " ";
    }
  }

  PredicatePtr boolean_predicate( KeywordPtr key ) {
    KeywordPredicate * pred = new KeywordPredicate(key);
    pred->addLiteral("1");
    return PredicatePtr(pred);
  }

  PredicatePtr same_keyword_predicate( KeywordPtr key, 
                                       const Selection& all,
                                       PredicatePtr pred ) {
    return PredicatePtr(new SamePredicate(key,all,pred));
  }

}}}

