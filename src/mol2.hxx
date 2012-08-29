#ifndef desres_msys_mol2_hxx
#define desres_msys_mol2_hxx

#include "system.hxx"

namespace desres { 
namespace msys {

    /* Write the given system as a mol2 file to the given stream */
    void ExportMol2( SystemPtr mol, std::ostream& out);

}}

#endif
