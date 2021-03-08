#ifndef desres_msys_analyze_hxx
#define desres_msys_analyze_hxx

#include "system.hxx"
#include <limits.h>
#include <chrono>

namespace desres { namespace msys {

    struct AssignBondOrder {
        enum Flags {
            Default = 0,
            ComputeResonantCharges = 1
        };
    };

    /* Assign bond order and formal charges to all fragments */
    void AssignBondOrderAndFormalCharge(SystemPtr mol, unsigned flags=0, std::chrono::milliseconds timeout=std::chrono::milliseconds(-1));

    /* Assign bond order and formal charges to the given atoms, all
     * of which should belong to the same fragment (i.e. they should
     * all be connected by bonds).  If total_charge is not supplied,
     * it will be guessed. */
    void AssignBondOrderAndFormalCharge(SystemPtr mol,
                                        IdList const& atoms,
                                        int total_charge = INT_MAX,
                                        unsigned flags = 0,
                                        std::chrono::milliseconds timeout=std::chrono::milliseconds(-1));

    /* Assign bond orders to aromatic bonds, leaving formal charges
     * and bonds between nonaromatic atoms alone.
     */
    void AssignAromaticBondOrders(SystemPtr mol);

    /* Compute topology ids.  FIXME definition please */
    IdList ComputeTopologicalIds(SystemPtr mol);

    /* Takes a molecule, and reorders the atoms so that all topologically
     * identical atoms are grouped together, and topological groups are sorted
     * by ascending id. The bonds are also sorted by the topological ids of their
     * constituent atoms.
     * returns a SystenPtr to the canonicalized system, and mappings
     * between the old atom/bond ids and the canonicalized ids
     */
    SystemPtr CanonicalizeMoleculeByTopids(SystemPtr mol,
                                           IdList const& atoms,
                                           std::map<Id, Id> &aid_to_canId,
                                           std::map<Id, Id> &bid_to_canId);

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

   /* Returns lists of non-pseudo bonds, pseudo bonds, angles, and dihedrals
    * in a given fragment or set of fragments */
    void GetBondsAnglesDihedrals(SystemPtr sys, const IdList& atoms,
            std::vector<IdList>& non_pseudo_bonds,
            std::vector<IdList>& pseudo_bonds,
            std::vector<IdList>& angles,
            std::vector<IdList>& dihedrals);


    /* check if the given set of atoms contains all its bonded neighbors.
     * The input set is assumed to be unique. */
    bool SelectionIsClosed(SystemPtr m, IdList const& ids, bool structure_only=false);
}}

#endif
