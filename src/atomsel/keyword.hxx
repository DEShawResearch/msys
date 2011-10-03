#ifndef desres_msys_atomsel_keyword_hxx
#define desres_msys_atomsel_keyword_hxx

#include "selection.hxx"
#include "regex.hxx"
#include <set>
#include <string>
#include <boost/shared_ptr.hpp>

namespace desres { namespace msys { namespace atomsel {

    /* an atomsel keyword returns values of a particular type.
     * There are three types: integer, floating point, and string,
     * represented by the C++ types int, double, and std::string.
     */
    enum KeywordType {
        KEY_INT=1, KEY_DBL=2, KEY_STR=4
    };
    typedef int Int;
    typedef double Dbl;
    typedef std::string Str;

    typedef std::string Literal;
    typedef std::pair<Literal,Literal> Range;

    struct Keyword {
        const std::string name;
        const KeywordType type;

        explicit Keyword( const std::string& n, KeywordType t ) 
            : name(n), type(t) {}
        virtual ~Keyword() {}

        /* get data of the specified type.  If the keyword isn't natively
         * of that type, the default implementations do the conversion.
         * Keyword subclasses MUST override the method corresponding to their
         * type. */
        virtual void iget( const Selection& s, std::vector<Int>& ) const;
        virtual void dget( const Selection& s, std::vector<Dbl>& ) const;
        virtual void sget( const Selection& s, std::vector<Str>& ) const;

        /* AND the literals and ranges into the given selection */
        void select( Selection& s,
                const std::set<Literal>& literals,
                const std::set<Range>&   ranges,
                const std::vector<Regex>& regexes) const;

    };

    typedef boost::shared_ptr<Keyword> KeywordPtr;

}}} // ns

#endif
