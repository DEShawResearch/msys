#ifndef desres_msys_sdf_hxx
#define desres_msys_sdf_hxx

#include "load.hxx"

namespace desres { namespace msys {

    /* Read the given MDL MOL/SDF file.  Return the first MOLECULE record */
    SystemPtr ImportSdf(std::string const& path);

    /* Iterator for SDF files */
    LoadIteratorPtr SdfIterator(std::string const& path);

    /* Read the given MDL MOL/SDF file, and return a list of Systems, one for
     * each MOLECULE record. */
    std::vector<SystemPtr> ImportSdfMany(std::string const& path);

    /* Write the structure to the given stream.  A single molecule
     * entry wil be created. */
    void ExportSdf( SystemPtr mol, std::ostream& out );

}}

#endif
