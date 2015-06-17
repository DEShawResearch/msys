#ifndef desres_msys_smiles_hxx
#define desres_msys_smiles_hxx

#include "molecule.hxx"

namespace desres { namespace msys {

    /* Construct a Molecule from a smiles string */
    MoleculePtr FromSmilesString(std::string const& smiles);

}}

#endif
