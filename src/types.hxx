#ifndef MOL_TYPES_HXX
#define MOL_TYPES_HXX

#include <string>
#include <set>
#include <vector>
#include <stdint.h>
#include <sstream>
#include <stdexcept>

#include "version.hxx"

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
    typedef std::vector<IdList> MultiIdList;

    enum { BadId = (uint32_t)-1 };
    inline bool bad(const Id& id) { return id==BadId; }

    struct Failure : public std::exception {
        explicit Failure(std::string const& msg) throw() : _msg(msg) {}
        virtual ~Failure() throw() {}
        virtual const char* what() const throw() { return _msg.c_str(); }

    private:
        std::string _msg;
    };

}}

#define MSYS_FAIL(args) do { \
    std::stringstream ss; \
    ss << args << "\nversion: " << MSYS_VERSION << "\nlocation: " << MSYS_LOC; \
    throw desres::msys::Failure(ss.str()); \
} while(0)

#endif
