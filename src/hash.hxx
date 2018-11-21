#ifndef desres_msys_hash_hxx
#define desres_msys_hash_hxx

#include "system.hxx"

namespace desres { namespace msys {

    // return a hash of the system, excluding provenance.
    uint64_t HashSystem(SystemPtr mol);

}}

#endif
