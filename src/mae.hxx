#ifndef desres_msys_mae_hxx
#define desres_msys_mae_hxx

#include "load.hxx"
#include <iostream>

namespace desres { namespace msys {

    LoadIteratorPtr MaeIterator(std::string const& path,
                                bool structure_only = false);

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized = false,
                         bool structure_only = false);

    SystemPtr ImportMAEFromStream( std::istream& in,
                         bool ignore_unrecognized = false,
                         bool structure_only = false);

    SystemPtr ImportMAEFromBytes( const char* bytes, int64_t len,
                         bool ignore_unrecognized = false,
                         bool structure_only = false);

    void ExportMAEMany( std::vector<SystemPtr> const& cts, 
                        std::string const& path,
                        bool with_forcefield = true );

    inline void ExportMAE( SystemPtr h, std::string const& path,
                    bool with_forcefield = true ) {
        ExportMAEMany(std::vector<SystemPtr>(1,h), path, with_forcefield);
    }
    
    
}}

#endif
