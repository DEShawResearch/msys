#include "ff.hxx"

using namespace desres::msys::mae;

namespace {
    typedef std::map<std::string, const Ffio *> FfioMap;
    FfioMap& registry() {
        static FfioMap r;
        return r;
    }
}

const Ffio * Ffio::get( const std::string& name ) {
    FfioMap::const_iterator i=registry().find(name);
    if (i!=registry().end()) return i->second;
    return NULL;
}

void Ffio::put( const std::string& name, const Ffio * f ) {
    registry()[name]=f;
}

