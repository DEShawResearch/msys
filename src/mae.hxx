#ifndef desres_msys_mae_hxx
#define desres_msys_mae_hxx

#include "system.hxx"

namespace desres { namespace msys {

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized = false );

    void ExportMAE( SystemPtr h, std::string const& path,
                    bool with_forcefield = true );
    
}}

#endif
