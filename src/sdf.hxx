#ifndef desres_msys_sdf_hxx
#define desres_msys_sdf_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Write the structure to the given stream.  A single molecule
     * entry wil be created. */
    void ExportSDF( SystemPtr mol, std::ostream& out );

}}

#endif
