#ifndef desres_msys_analyze_hxx
#define desres_msys_analyze_hxx

#include "system.hxx"
#include <limits.h>

namespace desres { namespace msys {

    struct AssignBondOrder {
        enum Flags {
            Default = 0,
            ComputeResonantCharges = 1
        };
    };

    /* Assign bond order and formal charges to all fragments */
    void AssignBondOrderAndFormalCharge(SystemPtr mol, unsigned flags=0);

    /* Assign bond order and formal charges to the given atoms, all
     * of which should belong to the same fragment (i.e. they should
     * all be connected by bonds).  If total_charge is not supplied,
     * it will be guessed. */
    void AssignBondOrderAndFormalCharge(SystemPtr mol,
                                        IdList const& atoms,
                                        int total_charge = INT_MAX,
                                        unsigned flags = 0);

    /* Assign bond orders to aromatic bonds, leaving formal charges
     * and bonds between nonaromatic atoms alone.
     */
    void AssignAromaticBondOrders(SystemPtr mol);

    /* Compute topology ids.  FIXME definition please */
    IdList ComputeTopologicalIds(SystemPtr mol);

    /* Add bonds based on estimated VDW radii determined from atomic
     * number (Bondi radii). */
    void GuessBondConnectivity(SystemPtr mol, bool periodic=false);

    /* Find representative fragments representing the complete set of
     * topologically distinct fragments, as determined by atomic number.
     *
     * Return mapping from representative fragment id to fragments having
     * identical topology.
     */
    std::map<Id,IdList> FindDistinctFragments(SystemPtr mol, MultiIdList const& fragments,
            std::vector<std::string> const& keys = std::vector<std::string>());

    /* Assign atom and residue types; do this after loading a new
     * system from a file or creating it from scratch.  This method
     * also calls updateFragids() for you. */
    void Analyze(SystemPtr mol);

    void GuessHydrogenPositions(SystemPtr mol, IdList const& hatoms);

}}

#endif
