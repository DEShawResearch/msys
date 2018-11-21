#ifndef desres_msys_analyze_get_fragments_hxx
#define desres_msys_analyze_get_fragments_hxx

#include "../system.hxx"

namespace desres { namespace msys {

    // Uses depth-first search to divide atomsIds list into connected components (fragments)
    void get_fragments(SystemPtr mol,
                       IdList const& atomIds, 
                       MultiIdList& fragments);

}}

#endif
