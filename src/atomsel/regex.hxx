#ifndef desres_msys_atomsel_regex_hxx
#define desres_msys_atomsel_regex_hxx

#include <string>
#include <boost/shared_ptr.hpp>

/* a wrapper class for pcre */

/* duplicate the forward declaration in pcre.h */
struct real_pcre;

namespace desres { namespace msys { namespace atomsel {

  class Regex {
  public:
    /* initialize with a regular expression pattern.  Throws on failure */
    Regex( const std::string& r, int pcre_options = 0 );

    /* does the input match? */
    bool match( const std::string& pat ) const;

  private:
    boost::shared_ptr<real_pcre> regex;
  };

}}} // ns

#endif
