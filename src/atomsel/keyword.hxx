#ifndef desres_msys_atomsel_keyword_hxx
#define desres_msys_atomsel_keyword_hxx

#include "selection.hxx"
#include "regex.hxx"
#include <set>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace desres { namespace msys { namespace atomsel {

    enum KeywordType {
        KEY_INT=1, KEY_DBL=2, KEY_STR=4
    };

    typedef int Int;
    typedef double Dbl;
    typedef std::string Str;

    struct Keyword : public boost::enable_shared_from_this<Keyword> {
        const KeywordType type;

        explicit Keyword(KeywordType t) 
        : type(t) {}
        virtual ~Keyword() {}

        virtual void iget(const Selection& s, std::vector<Int>& ) const;
        virtual void dget(const Selection& s, std::vector<Dbl>& ) const;
        virtual void sget(const Selection& s, std::vector<Str>& ) const;
    };

    typedef boost::shared_ptr<Keyword> KeywordPtr;

    struct TargetList {
        virtual ~TargetList() {}
        virtual void select(Selection& s, KeywordPtr k) = 0;
    };

    struct IntList : TargetList {
        typedef std::pair<Int,Int> Range;

        std::vector<Int>    values;
        std::vector<Range>  ranges;

        IntList() {}
        IntList(int i)          { add(i); }
        IntList(int i, int j)   { add(i,j); }
        void add(int i) { values.push_back(i); }
        void add(int i, int j) { ranges.push_back(Range(i,j)); }

        void select(Selection& s, KeywordPtr k);
    };

    struct FloatList : TargetList {
        typedef std::pair<Dbl,Dbl> Range;
        std::vector<Dbl>    values;
        std::vector<Range>  ranges;

        FloatList() {}
        FloatList(Dbl i) { add(i); }
        FloatList(Dbl i, Dbl j) { add(i,j); }
        void add(Dbl i) { values.push_back(i); }
        void add(Dbl i, Dbl j) { ranges.push_back(Range(i,j)); }

        void select(Selection& s, KeywordPtr k);
    };

    struct StringList : TargetList {
        std::vector<Str> values;
        std::vector<Regex*> regexes;

        StringList() {}
        StringList(Str const& s) { add(s); }
        StringList(Regex* r) { add(r); }
        void add(Str const& s) { values.push_back(s); }
        void add(Regex* r)  { regexes.push_back(r); }

        ~StringList() {
            for (unsigned i=0; i<regexes.size(); i++) delete regexes[i];
        }
        void select(Selection& s, KeywordPtr k);
    };

    Int parse_int(std::string const& s);
    Dbl parse_dbl(std::string const& s);

}}} // ns

#endif
