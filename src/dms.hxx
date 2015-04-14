#ifndef msys_dms_hxx
#define msys_dms_hxx

#include "system.hxx"
#include <iostream>

struct sqlite3;

namespace desres { namespace msys {

    // Changed in 1.7.101: structure_only includes pseudos, making it
    // consistent with export which also includes pseudos when
    // structure_only is true.
    SystemPtr ImportDMS(const std::string& path, bool structure_only=false);

    SystemPtr ImportDMSFromBytes( const char* bytes, int64_t len,
                                  bool structure_only=false);

    /* Return a hash of the contents of the dms file, or empty string if
     * the dms file does not contain hash information.  
     *
     * The hash excludes provenance and dms_version, in order to avoid
     * spurious change to the hash value resulting from simply reading the
     * file in and writing it back out again with no change.
     *
     * Unfortunately and unavoidably, the hash is sqlite version dependent.
     *
     * Note: a structure-only hash might be useful as well; however, some
     * thought should be given as to whether it ought to be a simply hash
     * of just the particle and bond tables, or else a hash of a 
     * canonicalized representation of the structure, like the InChI key.
     */
    String HashDMS(String const& path);

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

    namespace sqlite {
        SystemPtr ImportDMS(sqlite3* db, bool structure_only=false);
        void ExportDMS(SystemPtr sys, sqlite3* db,
                   Provenance const& provenance);
    }

}}

#endif
