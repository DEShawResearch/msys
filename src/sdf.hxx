#ifndef desres_msys_sdf_hxx
#define desres_msys_sdf_hxx

#include "io.hxx"

namespace desres { namespace msys {

    /* Read all entries from given SDF file */
    SystemPtr ImportSdf(std::string const& path);

    /* Iterator for SDF files */
    LoadIteratorPtr SdfFileIterator(std::string const& path);

    /* Iterator for SDF text */
    LoadIteratorPtr SdfTextIterator(std::string const& data);

    /* Export options */
    struct SdfExport {
        enum Flags {
            Default = 0,
            Append  = 1 << 0
        };
    };

    void ExportSdf(SystemPtr mol, std::string const& path, unsigned flags=0);
    std::string FormatSdf(SystemPtr mol);

}}

#endif
