#ifndef desres_msys_mol2_hxx
#define desres_msys_mol2_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Read the given mol2 file.  Return the first MOLECULE record */
    SystemPtr ImportMol2(std::string const& path);

    /* Read the given mol2 file, and return a list of Systems, one for
     * each MOLECULE record */
    std::vector<SystemPtr> ImportMol2Many(std::string const& path);

    /* Write the given system as a mol2 file to the given path.  A single
     * MOLECULE entry will be created.  Substructure ids will be numbered
     * from 1 based on the residue id.  No chain information will be 
     * written. */
    void ExportMol2( SystemPtr mol, std::string const& path,
                     Provenance const& provenance);

}}

#endif
