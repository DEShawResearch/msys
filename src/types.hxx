#ifndef MOL_TYPES_HXX
#define MOL_TYPES_HXX

#include <string>
#include <set>
#include <vector>
#include <stdint.h>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#ifdef _MSC_VER
#define MSYS_LOC __FILE__ << ":" << __LINE__ << "\n" << __FUNCSIG__
#else
#define MSYS_LOC __FILE__ << ":" << __LINE__ << "\n" << __PRETTY_FUNCTION__
#endif

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

    /* gettimeofday() */
    double now();

    const char* msys_version();
    int msys_major_version();
    int msys_minor_version();
    int msys_micro_version();
}}

#define MSYS_FAIL(args) do { \
    std::stringstream ss; \
    ss << args << "\nversion: " << desres::msys::msys_version() << "\nlocation: " << MSYS_LOC; \
    throw desres::msys::Failure(ss.str()); \
} while(0)

#endif
