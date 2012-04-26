#ifndef desres_msys_mae_hxx
#define desres_msys_mae_hxx

#include "system.hxx"
#include <iostream>

namespace desres { namespace msys {

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized = false );

    SystemPtr ImportMAEFromStream( std::istream& in,
                         bool ignore_unrecognized = false );

    SystemPtr ImportMAEFromBytes( const char* bytes, int64_t len,
                         bool ignore_unrecognized = false );

    void ExportMAE( SystemPtr h, std::string const& path,
                    bool with_forcefield = true );
    
}}

#endif
