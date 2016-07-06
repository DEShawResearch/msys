#ifndef desres_msys_smiles_hxx
#define desres_msys_smiles_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Construct a Molecule from a smiles string */
    SystemPtr FromSmilesString(std::string const& smiles,
                               bool forbid_stereo=true);

}}

#endif
