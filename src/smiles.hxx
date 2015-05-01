#ifndef desres_msys_smiles_hxx
#define desres_msys_smiles_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Construct a system from a smiles string */
    SystemPtr FromSmilesString(std::string const& smiles);

}}

#endif
