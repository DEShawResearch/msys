#ifndef desres_msys_sdf_hxx
#define desres_msys_sdf_hxx

#include "io.hxx"
#include "molecule.hxx"

namespace desres { namespace msys {

    /* Read the given MDL MOL/SDF file.  Return the first MOLECULE record */
    SystemPtr ImportSdf(std::string const& path);

    SystemPtr ImportSdfFromStream(std::istream& in);

    /* Iterator for SDF files */
    LoadIteratorPtr SdfIterator(std::string const& path);

    MoleculeIteratorPtr ScanSdf(FILE* fp);
    std::string FormatSdf( Molecule const& mol );

    struct SdfExport {
        enum Flags {
            Default = 0,
            Append  = 1 << 0
        };
    };

    /* Write the structure to the given stream.  A single molecule
     * entry wil be created. */
    void ExportSdf( SystemPtr mol, std::ostream& out );
    
    /* For commonality with other Export functions */
    void ExportSdf( SystemPtr mol, std::string const& path, unsigned flags=0);

}}

#endif
