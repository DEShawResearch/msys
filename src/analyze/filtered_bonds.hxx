#ifndef desres_viparr_filtered_bonds_hxx
#define desres_viparr_filtered_bonds_hxx

#include <msys/system.hxx>

namespace desres { namespace msys {

    /* Filter out virtuals/drudes to allow bond_order assigner to work
       with already typed systems */
    /* FIXME: Replace with msys>=1.5.12 provided alternatives */
    IdList filteredBondsForAtom(SystemPtr sys, Id aid);
    IdList filteredBondedAtoms(SystemPtr sys, Id aid);
    Id filteredBondCountForAtom(SystemPtr sys, Id aid);

}}

#endif
