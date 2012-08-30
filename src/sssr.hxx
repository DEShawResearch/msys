#ifndef desres_msys_sssr_hxx
#define desres_msys_sssr_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Find the smallest set of smallest rings for a fragment in a system.
     * The SSSR is not unique; if all_relevant is true, returns the union of 
     * all such sets. */
    MultiIdList GetSSSR(SystemPtr mol, IdList const& atoms,
            bool all_relevant=false);

    /* Helper functions and data structures for SSSR determination, exposed here
     * for testing/debugging purposes */
    namespace SSSR {

        /* We index vertices of a graph as 0,...,n-1. We represent an edge as an
         * index pair, a graph as a list of edges and a vertex-edge adjacency
         * list, and a subgraph as an indicator vector over edges plus an
         * optional vertex list (used to store an ordered path). */
        typedef std::pair<int, int> Edge;
        struct GraphRepr {
            std::vector<Edge> edges;
            std::vector<std::vector<int> > v_to_e;
            int other(int edge, int vertex) const {
                return ((edges[edge].first == vertex)
                        ? edges[edge].second : edges[edge].first);
            }
        };
        struct Subgraph {
            std::vector<bool> edges;
            std::vector<int> vertex_list;
            Subgraph() {}
            Subgraph(unsigned size, bool value) :
                edges(size, value), vertex_list() {}
            bool operator<(const Subgraph& other) const {
                return (vertex_list < other.vertex_list)
                    || ((vertex_list == other.vertex_list)
                            && (edges < other.edges));
            }
        };

        void get_biconnected_components(const GraphRepr& graph,
                std::vector<GraphRepr>& components,
                std::vector<std::vector<int> >& components_idx);
        void get_subgraph_path(const GraphRepr& graph,
                const Subgraph& subgraph, int start, int end, Subgraph& path,
                Subgraph& path_no_clear);
        void get_odd_path(const GraphRepr& graph, const Subgraph& odd_edges,
                int start, int end, bool return_multiple, int max_length,
                std::vector<Subgraph>& paths);
        unsigned get_cycle_basis(const GraphRepr& graph,
                std::deque<Subgraph>& basis, std::deque<int>& pivots,
                std::vector<int>& non_pivots);
        void minimize_cycle_basis(const GraphRepr& graph,
                std::deque<Subgraph>& basis, const std::deque<int>& pivots,
                unsigned nfixed, const std::vector<int>& non_pivots,
                std::vector<Subgraph>& min_basis);
        void get_relevant_cycles(const GraphRepr& graph,
                std::vector<Subgraph>& min_basis, const std::deque<int>& pivots,
                const std::vector<int>& non_pivots,
                std::set<Subgraph>& relevant_cycles);
    }

}}

#endif 
