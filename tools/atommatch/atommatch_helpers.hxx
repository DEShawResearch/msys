#ifndef desres_msys_atommatch_helpers_hxx
#define desres_msys_atommatch_helpers_hxx

#include <msys/system.hxx>
#include <msys/sssr.hxx>

namespace desres { namespace msys { namespace atommatch { namespace helpers {

    /* Create ScoreFct object that wraps C function pointer */
    class ScoreFctC : public ScoreFct {
        public:
            ScoreFctC(double (*c_apply)(SystemPtr, const IdList&, SystemPtr,
                        const IdList&)) : _c_apply(c_apply) { }
            virtual double apply(SystemPtr mol1, const IdList& atoms1,
                    SystemPtr mol2, const IdList& atoms2) const {
                return _c_apply(mol1, atoms1, mol2, atoms2);
            }
        private:
            double (*_c_apply)(SystemPtr, const IdList&, SystemPtr,
                    const IdList&);
    };

    /* Data types */
    using SSSR::Edge;
    using SSSR::GraphRepr;
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
    void init_graph(SystemPtr mol, std::vector<int>& atom_idx,
            GraphRepr& g);
    void partition_graph_biconnected(const GraphRepr& g,
            std::vector<GraphRepr>& components,
            std::vector<std::vector<int> >& components_idx,
            GraphRepr& tree, TreeEdges& tree_edges);
    void isomorphisms(const GraphRepr& g1,
            const std::vector<int>& comp_idx1,
            const std::vector<int>& atom_idx1, SystemPtr mol1,
            const GraphRepr& g2, const std::vector<int>& idx2,
            const std::vector<int>& atom_idx2, SystemPtr mol2,
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
            SystemPtr mol1,
            const std::vector<std::vector<int> >& components_idx2,
            const std::vector<int>& atom_idx2, SystemPtr mol2,
            const IsoMat& iso_mat, MatchList& match);

}}}}

#endif
