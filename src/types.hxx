#ifndef MOL_TYPES_HXX
#define MOL_TYPES_HXX

#include <string>
#include <set>
#include <vector>
#include <stdint.h>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream>

#ifdef _MSC_VER
#define MSYS_LOC __FILE__ << ":" << __LINE__ << "\n" << __FUNCSIG__
#else
#define MSYS_LOC __FILE__ << ":" << __LINE__ << "\n" << __PRETTY_FUNCTION__
#endif

#define MSYS_WARN(args) do { \
    std::cerr << args << "\nlocation: " << MSYS_LOC; \
} while(0)

#define MSYS_FAIL(args) do { \
    std::stringstream _msys_fail_tmp_ss_; \
    _msys_fail_tmp_ss_ << args << "\nlocation: " << MSYS_LOC; \
    throw desres::msys::Failure(_msys_fail_tmp_ss_.str()); \
} while(0)

namespace desres { namespace msys {

    typedef uint32_t Id;
    typedef int64_t  Int;
    typedef double   Float;
    typedef std::string String;
    typedef std::set<Id> IdSet;
    typedef std::vector<Id> IdList;
    typedef std::pair<Id,Id> IdPair;
    typedef std::vector<IdList> MultiIdList;

    static const Id BadId = -1;
    inline bool bad(const Id& id) { return id==BadId; }

    struct Failure : public std::exception {
        explicit Failure(std::string const& msg) throw() : _msg(msg) {}
        virtual ~Failure() throw() {}
        virtual const char* what() const throw() { return _msg.c_str(); }

    private:
        std::string _msg;
    };
    const char* msys_version();

    inline int stringToInt(std::string const& str) {
        char* stop;
        int res = strtol( str.c_str(), &stop, 10 );
        if ( *stop != 0 ) MSYS_FAIL("Bad int Specification: '" << str << "'");
        return res;
    }
    
    inline double stringToDouble(std::string const& str) {
        char* stop;
        double res = strtod( str.c_str(), &stop );
        if ( *stop != 0 ) MSYS_FAIL("Bad double Specification:\n" << str);
        return res;
    }


    static inline void trim(std::string& s) {
        size_t b=0, e=s.size();
        if (b==e) return;
        while (b<e && isspace(s[b])) ++b;
        --e;
        while (e>b && isspace(s[e])) --e;
        s=s.substr(b,e-b+1);
    }

    /* make the given container hold only unique elements in sorted order.
     * Return the number of non-unique elements (i.e. difference in size
     * between original and final container). */
    template <typename T>
    static Id sort_unique(T& t) {
        Id oldsize = t.size();
        std::sort(t.begin(), t.end());
        t.resize(std::unique(t.begin(), t.end()) - t.begin());
        return oldsize - t.size();
    }

    template <typename T>
    void to_lower(T& s) {
        std::for_each(std::begin(s), std::end(s), [](char& c) {c=tolower(c);});
    }
    template <typename T>
    void to_upper(T& s) {
        std::for_each(std::begin(s), std::end(s), [](char& c) {c=toupper(c);});
    }

    /* gettimeofday() */
    double now();

    int msys_major_version();
    int msys_minor_version();
    int msys_micro_version();
}}

#endif
