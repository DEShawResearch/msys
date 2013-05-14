#ifndef desres_msys_mae_hxx
#define desres_msys_mae_hxx

#include "io.hxx"
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

    struct MaeExport {
        enum Flags {
            Default             = 0,
            StructureOnly       = 1 << 0,
            CompressForcefield  = 1 << 1,
            Append              = 1 << 2
        };
    };

    void ExportMAEMany( std::vector<SystemPtr> const& cts, 
                        std::string const& path,
                        Provenance const& provenance,
                        unsigned flags=0);

    void ExportMAE( SystemPtr h, std::string const& path,
                           Provenance const& provenance,
                           unsigned flags=0);
    
}}

#endif
