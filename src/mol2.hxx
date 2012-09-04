#ifndef desres_msys_mol2_hxx
#define desres_msys_mol2_hxx

#include "system.hxx"

namespace desres { 
namespace msys {

    /* Write the given system as a mol2 file to the given path */
    void ExportMol2( SystemPtr mol, std::string const& path,
                     Provenance const& provenance);

}}

#endif
