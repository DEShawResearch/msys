#include "alchemical.hxx"
#include <sstream>
#include <stdexcept>
#include <boost/foreach.hpp>

using namespace desres::msys;

#define ERROR(x) do { \
    std::stringstream ss; \
    ss << "ERROR: " << x; \
    throw std::runtime_error(ss.str()); \
} while (0)



static void validate_map( IdList const& map, Id max, Id* nd ) {
    std::vector<Id> reals;
    Id d=0;
    BOOST_FOREACH(Id id, map) {
        if (bad(id)) {
            ++d; 
        } else if (id>=max) {
            ERROR("Invalid id " << id << " in atom map");
        } else {
            reals.push_back(id); 
        }
    }
    std::sort(reals.begin(), reals.end());
    reals.resize(std::unique(reals.begin(), reals.end())-reals.begin());
    if (reals.size() + d != map.size()) {
        ERROR("atom map contains duplicates, or some atoms not mapped");
    }
    if (reals.size() != max) {
        ERROR("not all atoms in system are mapped");
    }
    *nd = d;
}
                    

SystemPtr desres::msys::CreateAlchemical( SystemPtr A, IdList const& amap,
                                          SystemPtr B, IdList const& bmap ) {

    Id dumA, dumB;  /* number of unmapped atoms in amap and bmap */

    /* validate the atom map */
    try {
        validate_map(amap, A->maxAtomId(), &dumA);
    }
    catch (std::exception& e) {
        ERROR("Invalid map for A sites:\n" << e.what());
    }
    try {
        validate_map(bmap, B->maxAtomId(), &dumB);
    }
    catch (std::exception& e) {
        ERROR("Invalid map for B sites:\n" << e.what());
    }

    SystemPtr C;

    return C;
}
