#ifndef msys_dms_hxx
#define msys_dms_hxx

#include "system.hxx"

struct sqlite3;

namespace desres { namespace msys {

    SystemPtr ImportDMS(const std::string& path, bool structure_only=false);
    void ExportDMS(SystemPtr sys, const std::string& path, 
                   Provenance const& provenance);

    namespace sqlite {
        SystemPtr ImportDMS(sqlite3* db, bool structure_only=false);
        void ExportDMS(SystemPtr sys, sqlite3* db,
                   Provenance const& provenance);
    }

}}

#endif
