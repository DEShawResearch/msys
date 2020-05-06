#ifndef desres_msys_psf_hxx
#define desres_msys_psf_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr ImportPSF( std::string const& path );
    void ExportPSF(SystemPtr mol, std::string const& path);

}}

#endif
