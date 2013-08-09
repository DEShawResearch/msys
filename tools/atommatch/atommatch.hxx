#ifndef desres_msys_atommatch_hxx
#define desres_msys_atommatch_hxx

#include <msys/system.hxx>
#include <boost/shared_ptr.hpp>

namespace desres { namespace msys { namespace atommatch {

    /* Interface for a function that assigns numeric scores to isomorphisms
     * between biconnected components */
    struct ScoreFct {
        virtual double apply(SystemPtr mol1, const IdList& atoms1,
                SystemPtr mol2, const IdList& atoms2) const = 0;
        virtual ~ScoreFct() { }
    };
    typedef boost::shared_ptr<ScoreFct> ScoreFctPtr;

    /* AtomMatch
     *
     * Arguments:
     *   mol1 -- first system with a single molecule
     *   mol2 -- second system with a single molecule
     *   score_fct -- user-defined function to score component isomorphisms
     *   atom1 -- optional atom in mol1 that must be matched to atom2
     *   atom2 -- optional atom in mol2 that must be matched to atom1
     *
     * Returns:
     *   match -- list of ID pairs where the first ID is an atom in mol1 and
     *       the second is the matching atom in mol2
     *
     * This function returns the "best" isomorphism between subgraphs in mol1
     * and mol2, where "best" is determined by the following rules:
     * (1) Each subgraph is a union of complete biconnected components of the
     *     molecule (so we cannot match only part of a ring system).
     * (2) The match has the highest score among all possible subgraph matches,
     *     where the score of a match is the sum of scores over matching
     *     biconnected components as determined by 'score_fct'.
     * (3) Among multiple matches satisfying (2), the match contains the most
     *     matched atoms.
     * (4) Among multiple matches satisfying (3), the match minimizes RMSD
     *     between the matched atoms.
     *
     * Rules (1)-(3) are enforced exactly. As the number of matches that satisfy
     * rules (1)-(3) may increase exponentially in the number of
     * symmetries of the molecule, rule (4) is not enforced exactly but rather
     * is implemented by a heuristic search.
     */
    typedef std::vector<std::pair<int, int> > MatchList;
    void AtomMatch(SystemPtr mol1, SystemPtr mol2, ScoreFctPtr score_fct,
            MatchList& match, Id atom1=BadId, Id atom2=BadId);

}}}

#endif
