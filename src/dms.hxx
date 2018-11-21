#ifndef msys_dms_hxx
#define msys_dms_hxx

#include "system.hxx"
#include "io.hxx"
#include <iostream>

struct sqlite3;

namespace desres { namespace msys {

    LoadIteratorPtr DMSIterator(std::string const& path,
                                bool structure_only,
                                bool without_tables);

    // Changed in 1.7.101: structure_only includes pseudos, making it
    // consistent with export which also includes pseudos when
    // structure_only is true.
    //
    // Changed in 1.7.109: revert this change and introduce without_tables.
    SystemPtr ImportDMS(const std::string& path, bool structure_only,
                                                 bool without_tables);

    // legacy API
    inline
    SystemPtr ImportDMS(const std::string& path, bool structure_only=false) {
        return ImportDMS(path, structure_only, structure_only);
    }


    SystemPtr ImportDMSFromBytes( const char* bytes, int64_t len,
                                  bool structure_only, bool without_tables);

    inline
    SystemPtr ImportDMSFromBytes( const char* bytes, int64_t len,
                                  bool structure_only=false) {
        return ImportDMSFromBytes(bytes, len, structure_only, structure_only);
    }


    struct DMSExport {
        enum Flags { Default            = 0 
                   , Append             = 1 << 0
                   , StructureOnly      = 1 << 1
                   , Unbuffered         = 1 << 2
        };
    };

    void ExportDMS(SystemPtr sys, const std::string& path, 
                   Provenance const& provenance,
                   unsigned flags = 0);

    std::string FormatDMS(SystemPtr sys, Provenance const& prov);

    namespace sqlite {
        SystemPtr ImportDMS(sqlite3* db, bool structure_only, bool with_tables);

        inline
        SystemPtr ImportDMS(sqlite3* db, bool structure_only=false) {
            return ImportDMS(db, structure_only, structure_only);
        }
        void ExportDMS(SystemPtr sys, sqlite3* db,
                   Provenance const& provenance);
    }

}}

#endif
