#ifndef desres_msys_pdb_hxx
#define desres_msys_pdb_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr ImportPDB( std::string const& path );
    void ExportPDB(SystemPtr mol, std::string const& path);

}}

#endif
