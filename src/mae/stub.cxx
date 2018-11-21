#include "../mae.hxx"

namespace desres { namespace msys {

    SystemPtr ImportMAE( std::string const& path,
                         bool ignore_unrecognized ) {
        MSYS_FAIL("no mae file support");
    }

    SystemPtr ImportMAEFromStream( std::istream& in,
                         bool ignore_unrecognized ) {
        MSYS_FAIL("no mae file support");
    }

    SystemPtr ImportMAEFromBytes( const char* bytes, int64_t len,
                         bool ignore_unrecognized ) {
        MSYS_FAIL("no mae file support");
    }

    void ExportMAE( SystemPtr h, std::string const& path,
                    bool with_forcefield ) {
        MSYS_FAIL("no mae file support");
    }
    
}}

