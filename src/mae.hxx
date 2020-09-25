#ifndef desres_msys_mae_hxx
#define desres_msys_mae_hxx

#include "io.hxx"
#include <iostream>

namespace desres { namespace msys {

    LoadIteratorPtr MaeIterator(std::string const& path,
                                bool structure_only = false);

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized,
                         bool structure_only,
                         bool without_tables);

    inline
    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized = false,
                         bool structure_only = false) {
        return ImportMAE(path, ignore_unrecognized, 
                                structure_only, structure_only);
    }

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
            //CompressForcefield  = 1 << 1, /* deprecated */
            Append              = 1 << 2,
            AllowReorderAtoms   = 1 << 3,
        };
    };

    std::string ExportMAEContents( SystemPtr h,
                            Provenance const& provenance,
                            unsigned flags=0);

    void ExportMAE( SystemPtr h, std::string const& path,
                           Provenance const& provenance,
                           unsigned flags=0);

}}

#endif
