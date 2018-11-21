#ifndef desres_msys_xyz_hxx
#define desres_msys_xyz_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr ImportXYZ( std::string const& path );
    void ExportXYZ( SystemPtr mol, std::string const& path );

}}

#endif
