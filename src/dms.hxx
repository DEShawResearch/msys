#ifndef msys_dms_hxx
#define msys_dms_hxx

#include "system.hxx"
#include <iostream>

struct sqlite3;

namespace desres { namespace msys {

    SystemPtr ImportDMS(const std::string& path, bool structure_only=false);

    SystemPtr ImportDMSFromBytes( const char* bytes, int64_t len,
                                  bool structure_only=false);

    struct DMSExport {
        enum Flags { Default            = 0 
                   , Append             = 1 << 0
                   , StructureOnly      = 1 << 1
        };
    };

    void ExportDMS(SystemPtr sys, const std::string& path, 
                   Provenance const& provenance,
                   unsigned flags = 0);

    namespace sqlite {
        SystemPtr ImportDMS(sqlite3* db, bool structure_only=false);
        void ExportDMS(SystemPtr sys, sqlite3* db,
                   Provenance const& provenance);
    }

}}

#endif
