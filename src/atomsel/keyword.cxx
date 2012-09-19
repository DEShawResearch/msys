#include "keyword.hxx"
#include <cstdlib>
#include <stdexcept>
#include <errno.h>
#include <boost/static_assert.hpp>

using namespace desres::msys::atomsel;

Int desres::msys::atomsel::parse_int(std::string const& s) {
    static const int base=10;
    char* end;
    errno=0;
    Int result = strtoll(s.c_str(), &end, base);
    if (errno==ERANGE) {
        MSYS_FAIL("Integer value " << s << " would overflow");
    }
    if (s.empty() || *end!='\0') {
        MSYS_FAIL("Could not convert all of '" << s << "' to integer");
    }
    return result;
}

Dbl desres::msys::atomsel::parse_dbl(std::string const& s) {
    char* end;
    errno=0;
    Dbl result = strtod(s.c_str(), &end);
    if (errno==ERANGE) {
        MSYS_FAIL("Float value " << s << " would overflow");
    }
    if (s.empty() || *end!='\0') {
        MSYS_FAIL("Could not convert all of '" << s << "' to float");
    }
    return result;
}

void Keyword::iget( const Selection& s, std::vector<Int>& v ) const {
  int i,n=s.size();
  switch (type) {
    case KEY_INT:
      throw std::runtime_error("Expected override of iget");
      break;
    case KEY_DBL:
      {
        std::vector<Dbl> tmp(n);
        dget(s,tmp);
        for (i=0; i<n; i++) v[i]=static_cast<Int>(tmp[i]);
      }
      break;
    case KEY_STR:
      for (i=0; i<n; i++) v[i]=0;
      break;
    default:
      throw std::runtime_error("Unsupported key type");
  }
}

void Keyword::dget( const Selection& s, std::vector<Dbl>& v ) const {
  int i,n=s.size();
  switch (type) {
    case KEY_INT:
      {
        std::vector<Int> tmp(n);
        iget(s,tmp);
        for (i=0; i<n; i++) v[i]=tmp[i];
      }
      break;
    case KEY_DBL:
      throw std::runtime_error("Expected override of dget");
      break;
    case KEY_STR:
      for (i=0; i<n; i++) v[i]=0;
      break;
    default:
      throw std::runtime_error("Unsupported key type");
  }
}

void Keyword::sget( const Selection& s, std::vector<Str>& v ) const {
  int i,n=s.size();
  switch (type) {
    case KEY_INT:
    case KEY_DBL:
      for (i=0; i<n; i++) v[i]="";
      break;
    case KEY_STR:
      throw std::runtime_error("Expected override of sget");
      break;
    default:
      throw std::runtime_error("Unsupported key type");
  }
}

namespace {
  template <typename T>
    void compare_literals( const std::vector<T>& v,
        const std::set<T>& lits,
        Selection& s ) {
      int i,n=s.size();
      for (i=0; i<n; i++) s[i] |= lits.count(v[i]);
    }

  template <typename T>
    void compare_ranges( const std::vector<T>& v,
        const std::set<std::pair<T,T> >& ranges,
        Selection& s ) {
      int i,n=s.size();
      for (i=0; i<n; i++) {
        if (s[i]) continue;
        typedef typename std::set<std::pair<T,T> >::const_iterator iter;
        for (iter r=ranges.begin(); r!=ranges.end(); ++r) {
          const T& min=r->first;
          const T& max=r->second;
          s[i] |= (v[i]>=min && v[i]<=max);
        }
      }
    }

  void compare_regexes( const std::vector<std::string>& v,
                        const std::vector<Regex>& regexes,
                        Selection& s ) {
    int i,n=s.size();
    for (i=0; i<n; i++) {
      if (s[i]) continue;
      for (unsigned j=0; j<regexes.size(); j++) {
        s[i] |= regexes[j].match(v[i]);
      }
    }
  }
}

void Keyword::select( Selection& s,
    const std::set<Literal>& literals,
    const std::set<Range>&   ranges,
    const std::vector<Regex>& regexes) const {

  Selection s2(s);
  s2.clear();
  switch (type) {
    case KEY_INT:
      {
        if (regexes.size())
            MSYS_FAIL("cannot select on regex for keyword '" << name << "' of integer type");

        std::vector<Int> v(s.size());
        iget(s,v);
        std::set<Int> lit;
        std::set<std::pair<Int,Int> > ran;
        for (std::set<Literal>::const_iterator i=literals.begin();
            i!=literals.end(); i++) {
          lit.insert(parse_int(*i));
        }
        for (std::set<Range>::const_iterator i=ranges.begin();
            i!=ranges.end(); i++) {
          ran.insert(std::make_pair(parse_int(i->first),parse_int(i->second)));
        }
        compare_literals( v, lit, s2 );
        compare_ranges(   v, ran, s2 );
      }
      break;
    case KEY_DBL:
      {
        if (regexes.size())
          MSYS_FAIL("cannot select on regex for keyword '" << name << "' of float type");

        std::vector<Dbl> v(s.size());
        dget(s,v);
        std::set<Dbl> lit;
        std::set<std::pair<Dbl,Dbl> > ran;
        for (std::set<Literal>::const_iterator i=literals.begin();
            i!=literals.end(); i++) {
          lit.insert(parse_dbl(*i));
        }
        for (std::set<Range>::const_iterator i=ranges.begin();
            i!=ranges.end(); i++) {
          ran.insert(std::make_pair(parse_dbl(i->first),parse_dbl(i->second)));
        }
        compare_literals( v, lit, s2 );
        compare_ranges(   v, ran, s2 );
      }
      break;
    case KEY_STR:
      {
        std::vector<Str> v(s.size());
        sget(s,v);
        compare_literals( v, literals, s2 );
        compare_ranges(   v, ranges,   s2 );
        compare_regexes(  v, regexes,  s2 );
      }
      break;
    default:
      throw std::runtime_error("Unsupported key type");
  }
  s.intersect(s2);
}
