#ifndef MOL_TYPES_HXX
#define MOL_TYPES_HXX

#include <string>
#include <set>
#include <vector>
#include <stdint.h>

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

}}

#endif
