#ifndef desres_msys_analyze_hxx
#define desres_msys_analyze_hxx

#include "system.hxx"
#include <limits.h>

namespace desres { namespace msys {

    /* Assign bond order and formal charges to all fragments */
    void AssignBondOrderAndFormalCharge(SystemPtr mol);

    /* Assign bond order and formal charges to the given atoms, all
     * of which should belong to the same fragment (i.e. they should
     * all be connected by bonds).  If total_charge is not supplied,
     * it will be guessed. */
    void AssignBondOrderAndFormalCharge(SystemPtr mol,
                                        IdList const& atoms,
                                        int total_charge = INT_MAX);

    /* Compute topology ids.  FIXME definition please */
    IdList ComputeTopologicalIds(SystemPtr mol);

    /* Add bonds based on estimated VDW radii determined from atomic
     * number (Bondi radii). */
    void GuessBondConnectivity(SystemPtr mol);

    /* Find representative fragments representing the complete set of
     * topologically distinct fragments, as determined by atomic number. */
    IdList FindDistinctFragments(SystemPtr mol, MultiIdList const& fragments);

    void GuessHydrogenPositions(SystemPtr mol, IdList const& hatoms);
}}

#endif
