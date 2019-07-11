#ifndef desres_msys_json_hxx
#define desres_msys_json_hxx

#include "io.hxx"

namespace desres { namespace msys {

    SystemPtr ImportJson(std::string const& path);
    SystemPtr ParseJson(const char* text);

    struct JsonExport {
        enum Flags {
              Default       = 0
            , StructureOnly = 1 << 1
            , Whitespace    = 1 << 2
        };
    };

    void ExportJson(SystemPtr mol, std::string const& path, Provenance const& provenance,
            unsigned flags=0);

    std::string FormatJson(SystemPtr mol, Provenance const& provenance, unsigned flags=0);

}}

#endif
