#ifndef desres_msys_mol2_hxx
#define desres_msys_mol2_hxx

#include "io.hxx"

namespace desres { namespace msys {

    /* Read the given mol2 file.  Return the first MOLECULE record */
    SystemPtr ImportMol2(std::string const& path);

    /* Read the given mol2 file, and return a list of Systems, one for
     * each MOLECULE record.  The atom type will be stored in the atom
     * property "sybyl_type", and the bond type in the bond property
     * "sybyl_type".  */
    std::vector<SystemPtr> ImportMol2Many(std::string const& path);

    /* Iterator for mol2 files */
    LoadIteratorPtr Mol2Iterator(std::string const& path);

    /* Assign sybyl atom and bond types to the given system.  Be sure
     * bond order is valid; see AssignBondOrderAndFormalCharge.
     */
    const char* GuessSybylAtomType(SystemPtr mol, Id id, bool cyclic);
    const char* GuessSybylBondType(SystemPtr mol, Id bond_id, 
                                   std::string const& itype,
                                   std::string const& jtype);

    struct Mol2Export {
        enum Flags {
            Default     = 0,
            Append      = 1 << 0
        };
    };

    /* Write the given system as a mol2 file to the given path.  A single
     * MOLECULE entry will be created.  Substructure ids will be numbered
     * from 1 based on the residue id.  No chain information will be 
     * written.  If "sybyl_type" is found in atom properties, it will be
     * used for the atom type; otherwise, the atom type will be identical
     * to the chemical element abbreviation.  If "sybyl_type" is found in
     * bond properties, it will be used for the bond type; otherwise,
     * the bond type will be "1", "2", etc. according to the bond order. */
    void ExportMol2( SystemPtr mol, std::string const& path,
                     Provenance const& provenance,
                     IdList const& ids,
                     unsigned flags=0);

}}

#endif
