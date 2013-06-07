#include <msys/system.hxx>
#include <msys/sssr.hxx>
#include <boost/shared_ptr.hpp>
#include <stdexcept>

#define FAIL(str) throw std::runtime_error(str)

namespace desres { namespace fep_atommatch {

    /* Interface for a function that assigns numeric scores to isomorphisms
     * between biconnected components */
    struct ScoreFct {
        virtual double apply(msys::SystemPtr mol1, const msys::IdList& atoms1,
                msys::SystemPtr mol2, const msys::IdList& atoms2) const = 0;
        virtual ~ScoreFct() { }
    };
    typedef boost::shared_ptr<ScoreFct> ScoreFctPtr;

    /* FEPAtomMatch
     *
     * Arguments:
     *   mol1 -- first system with a single molecule
     *   mol2 -- second system with a single molecule
     *   score_fct -- user-defined function to score component isomorphisms
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
    void FEPAtomMatch(msys::SystemPtr mol1, msys::SystemPtr mol2,
            ScoreFctPtr score_fct, MatchList& match);

    /* Helper functions and data structures, exposed here for testing/debugging
     * purposes */
    namespace helpers {
        /* Create ScoreFct object that wraps C function pointer */
        class ScoreFctC : public ScoreFct {
            public:
                ScoreFctC(double (*c_apply)(msys::SystemPtr,
                            const msys::IdList&, msys::SystemPtr,
                            const msys::IdList&))
                    : _c_apply(c_apply) { }
                virtual double apply(msys::SystemPtr mol1,
                        const msys::IdList& atoms1,
                        msys::SystemPtr mol2,
                        const msys::IdList& atoms2) const {
                    return _c_apply(mol1, atoms1, mol2, atoms2);
                }
            private:
                double (*_c_apply)(msys::SystemPtr, const msys::IdList&,
                        msys::SystemPtr, const msys::IdList&);
        };

        /* Data types */
        using msys::SSSR::Edge;
        using msys::SSSR::GraphRepr;
        struct Isomorphism {
            std::vector<int> perm;
            double score;
        };
        typedef std::vector<Isomorphism> Isomorphisms;
        typedef std::vector<std::vector<Isomorphisms> > IsoMat;
        typedef std::vector<std::vector<std::vector<std::pair<int, int> > > >
            TreeEdges;
        struct DPNode {
            int v1;
            int v2;
            int iso_ind;
            DPNode(int _v1, int _v2, int _iso_ind)
                : v1(_v1), v2(_v2), iso_ind(_iso_ind) { }
        };
        typedef std::vector<DPNode> DPNodes;
        struct DPBest {
            std::vector<DPNodes> options;
            double score;
            int size;
        };
        typedef std::vector<std::vector<std::vector<DPBest> > > DPMat;
 
        /* Helper functions, documented in cxx file */
        void init_graph(msys::SystemPtr mol, std::vector<int>& atom_idx,
                GraphRepr& g);
        void partition_graph_biconnected(const GraphRepr& g,
                std::vector<GraphRepr>& components,
                std::vector<std::vector<int> >& components_idx,
                GraphRepr& tree, TreeEdges& tree_edges);
        void isomorphisms(const GraphRepr& g1,
                const std::vector<int>& comp_idx1,
                const std::vector<int>& atom_idx1, msys::SystemPtr mol1,
                const GraphRepr& g2, const std::vector<int>& idx2,
                const std::vector<int>& atom_idx2, msys::SystemPtr mol2,
                ScoreFctPtr rep, Isomorphisms& isos);
        void get_levels(const GraphRepr& tree, int root,
                std::vector<std::vector<int> >& levels, std::vector<int>&
                level_map, bool screen);
        void get_perms(int size1, int size2,
                std::vector<std::vector<int> >& perms);
        void options_append(std::vector<DPNodes>& options, int v1, int v2,
                const std::vector<int>& option_inds);
        void dp_update(int v1, int v2, DPMat& best_matches,
                const std::vector<std::vector<std::pair<int, int> > >&
                tree_edges1, const std::vector<int>& level_map1,
                const std::vector<std::vector<std::pair<int, int> > >&
                tree_edges2, const std::vector<int>& level_map2,
                const IsoMat& iso_mat, bool keep_options=true);
        std::pair<double, int> max_match(const GraphRepr& tree1,
                const TreeEdges& tree_edges1,
                const std::vector<GraphRepr>& components1,
                const std::vector<std::vector<int> > components_idx1, int root1,
                const GraphRepr& tree2, const TreeEdges& tree_edges2,
                const std::vector<GraphRepr>& components2,
                const std::vector<std::vector<int> > components_idx2, int root2,
                const IsoMat& iso_mat, DPMat& best_matches,
                bool keep_options=true);
        double select_single(const DPMat& best_matches, int v1, int v2,
                int iso_ind, const std::vector<std::vector<int> >&
                components_idx1, const std::vector<int>& atom_idx1,
                msys::SystemPtr mol1,
                const std::vector<std::vector<int> >& components_idx2,
                const std::vector<int>& atom_idx2, msys::SystemPtr mol2,
                const IsoMat& iso_mat, MatchList& match);
    }

}}
