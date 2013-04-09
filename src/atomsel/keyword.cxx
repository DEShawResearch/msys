#include "keyword.hxx"
#include <cstdlib>
#include <stdexcept>
#include <errno.h>

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
                           const std::vector<T>& lits,
                           Selection& s ) {
      int i,n=s.size();
      for (i=0; i<n; i++) {
          s[i] |= std::binary_search(lits.begin(), lits.end(), v[i]);
      }
    }

  template <typename T>
    void compare_ranges( const std::vector<T>& v,
                         const std::vector<std::pair<T,T> >& ranges,
                           Selection& s ) {
      int i,n=s.size();
      for (i=0; i<n; i++) {
        if (s[i]) continue;
        typedef typename std::vector<std::pair<T,T> >::const_iterator iter;
        for (iter r=ranges.begin(); r!=ranges.end(); ++r) {
          const T& min=r->first;
          const T& max=r->second;
          s[i] |= (v[i]>=min && v[i]<=max);
        }
      }
    }

  void compare_regexes( const std::vector<std::string>& v,
                        const std::vector<Regex*>& regexes,
                        Selection& s ) {
    int i,n=s.size();
    for (i=0; i<n; i++) {
      if (s[i]) continue;
      for (unsigned j=0; j<regexes.size(); j++) {
        s[i] |= boost::regex_match(v[i], *regexes[j]);
      }
    }
  }
}

void IntList::select(Selection& s, KeywordPtr k) {
    /* fetch values for current selection */
    std::vector<Int> v(s.size());
    k->iget(s,v);

    /* compute the size of a perfect hash table */
    Int imin=0, imax=0;
    bool have_first = false;
    for (Id i=0; i<values.size(); i++) {
        Int const& v = values[i];
        if (have_first) {
            imin = std::min(imin,v);
            imax = std::max(imax,v);
        } else {
            have_first = true;
            imin = v;
            imax = v;
        }
    }

    /* add selected values */
    std::vector<bool> table(1+imax-imin);
    for (Id i=0; i<values.size(); i++) {
        Int const& v = values[i];
        if (v >= imin && v<= imax) {
            table[v-imin] = true;
        }
    }

    /* update selection */
    for (Id i=0; i<s.size(); i++) {
        if (s[i]) {
            Int x = v[i];
            s[i] = x >= imin && x<= imax && table[x-imin];
            if (!s[i]) {
                bool on = false;
                for (Id j=0; j<ranges.size(); j++) {
                    if (x >= ranges[j].first && x <= ranges[j].second) {
                        on = true;
                        break;
                    }
                }
                s[i] = on;
            }
        }
    }
}

void FloatList::select(Selection& s, KeywordPtr k) {
    /* FIXME lots of extra work here */
    Selection s2(s);
    s2.clear();
    std::vector<Dbl> v(s.size());
    k->dget(s,v);
    sort_unique(values);
    compare_literals(v, values, s2);
    compare_ranges(  v, ranges, s2);
    s.intersect(s2);
}

void StringList::select(Selection& s, KeywordPtr k) {
    /* FIXME lots of extra work here */
    Selection s2(s);
    s2.clear();
    std::vector<Str> v(s.size());
    k->sget(s,v);
    sort_unique(values);
    compare_literals(v, values, s2);
    //compare_ranges(  v, ran, s2);
    compare_regexes( v, regexes,s2 );
    s.intersect(s2);
}

