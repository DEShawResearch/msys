#ifndef desres_msys_mae_hxx
#define desres_msys_mae_hxx

#include "system.hxx"
#include <iostream>

namespace desres { namespace msys {

    std::vector<SystemPtr> ImportMAEMany(std::string const& path);
    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized = false,
                         bool structure_only = false);

    SystemPtr ImportMAEFromStream( std::istream& in,
                         bool ignore_unrecognized = false,
                         bool structure_only = false);

    SystemPtr ImportMAEFromBytes( const char* bytes, int64_t len,
                         bool ignore_unrecognized = false,
                         bool structure_only = false);

    void ExportMAE( SystemPtr h, std::string const& path,
                    bool with_forcefield = true );
    
}}

#endif
