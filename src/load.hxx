#ifndef msys_load_hxx
#define msys_load_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Try to guess the file type from the path, and load the system
     * using default options */
    SystemPtr Load(std::string const& path);

}}

#endif
