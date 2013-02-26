#include "../sssr.hxx"
#include "bondFilters.hxx"
#include <stack>
#include <queue>
#include <set>
#include <algorithm>

using namespace desres::msys::SSSR;

namespace {
    /* Hopcroft-Tarjan algorithm for biconnected components. See R. Tarjan,
     * "Depth-first search and linear graph algorithms", SIAM J. Comput. 1(2),
     * 1972. */
    struct BCCItem { /* Vertex in DFS */
        int idx; /* Vertex index */
        int lowpoint; /* Depth of shallowest neighbor of some descendant */
        int parent; /* Parent in DFS */
        int next_nbr; /* Next neighbor to check in DFS */
        BCCItem(int _idx, int _lowpoint, int _parent) :
            idx(_idx), lowpoint(_lowpoint), parent(_parent), next_nbr(0) {}
    };
}

void desres::msys::SSSR::get_biconnected_components(const GraphRepr& graph,
        std::vector<GraphRepr>& components,
        std::vector<std::vector<int> >& components_idx) {

    components.clear();
    components_idx.clear();
    
    /* Depth of vertices in DFS */
    std::vector<int> depths(graph.v_to_e.size(), -1);

    /* Loop over connected components */
    for (unsigned root = 0; root < graph.v_to_e.size(); ++root) {
        /* Ignore empty vertex indices and vertices that have already been
         * processed */
        if (graph.v_to_e[root].size() == 0) continue;
        if (depths[root] >= 0) continue;
        /* Stack of vertices (and other info in BCCItem) maintained in a DFS */
        std::stack<BCCItem> DFS_stack;
        /* Stack of directed edges, represented by a
         * (tail_vertex_idx, edge_idx) pair */
        std::stack<std::pair<int, int> > edge_stack;

        /* Do DFS */
        DFS_stack.push(BCCItem(root,0,-1));
        depths[root] = 0;
        while (DFS_stack.size() > 0) {
            BCCItem& item = DFS_stack.top();
            const std::vector<int>& edges = graph.v_to_e[item.idx];
            bool found = false;
            for (unsigned j = item.next_nbr; j < edges.size(); ++j) {
                ++item.next_nbr;
                int nbr = graph.other(edges[j], item.idx);
                if (depths[nbr] == -1) {
                    /* Found unexplored child in DFS; add tree edge and add
                     * child to DFS tree */
                    edge_stack.push(std::make_pair(item.idx, edges[j]));
                    DFS_stack.push(BCCItem(nbr, depths[item.idx]+1, item.idx));
                    depths[nbr] = depths[item.idx]+1;
                    found = true;
                    break;
                } else if (depths[nbr] < depths[item.idx]
                        && nbr != item.parent) {
                    /* Neighbor is an ancestor; add backwards edge and update
                     * lowpoint */
                    edge_stack.push(std::make_pair(item.idx, edges[j]));
                    item.lowpoint = (item.lowpoint < depths[nbr] 
                            ? item.lowpoint : depths[nbr]);
                }
            }

            /* If we found an unexplored child, continue with that child in the
             * DFS */
            if (found) continue;

            int lowpoint = item.lowpoint;
            /* Remove current vertex from DFS stack */
            DFS_stack.pop();
            if (DFS_stack.size() > 0) {
                /* Update low point of current parent */
                DFS_stack.top().lowpoint = (DFS_stack.top().lowpoint < lowpoint
                        ? DFS_stack.top().lowpoint : lowpoint);
                if (lowpoint >= depths[DFS_stack.top().idx]) {
                    /* Found biconnected component; create new GraphRepr for
                     * component and a 2-way vertex index mapping */
                    GraphRepr comp;
                    std::vector<int> comp_to_graph;
                    std::map<int, int> graph_to_comp;

                    /* Pop off edges and add to component while the head vertex
                     * is deeper than the current parent */
                    while (edge_stack.size() > 0
                            && depths[edge_stack.top().first] >
                            depths[DFS_stack.top().idx]) {
                        int v1 = edge_stack.top().first;
                        std::map<int, int>::iterator iter1
                            = graph_to_comp.find(v1);
                        if (iter1 == graph_to_comp.end()) {
                            comp.v_to_e.push_back(std::vector<int>());
                            comp_to_graph.push_back(v1);
                            iter1 = graph_to_comp.insert(std::make_pair(v1,
                                        comp_to_graph.size()-1)).first;
                        }
                        int v2 = graph.other(edge_stack.top().second, v1);
                        std::map<int, int>::iterator iter2
                            = graph_to_comp.find(v2);
                        if (iter2 == graph_to_comp.end()) {
                            comp.v_to_e.push_back(std::vector<int>());
                            comp_to_graph.push_back(v2);
                            iter2 = graph_to_comp.insert(std::make_pair(v2,
                                        comp_to_graph.size()-1)).first;
                        }
                        comp.edges.push_back(Edge(iter1->second,iter2->second));
                        comp.v_to_e[iter1->second].push_back(
                                comp.edges.size()-1);
                        comp.v_to_e[iter2->second].push_back(
                                comp.edges.size()-1);
                        edge_stack.pop();
                    }
                    /* Pop off edge and add to component one more time */
                    int v1 = edge_stack.top().first;
                    std::map<int, int>::iterator iter1 = graph_to_comp.find(v1);
                    if (iter1 == graph_to_comp.end()) {
                        comp.v_to_e.push_back(std::vector<int>());
                        comp_to_graph.push_back(v1);
                        iter1 = graph_to_comp.insert(std::make_pair(v1,
                                    comp_to_graph.size()-1)).first;
                    }
                    int v2 = graph.other(edge_stack.top().second, v1);
                    std::map<int, int>::iterator iter2 = graph_to_comp.find(v2);
                    if (iter2 == graph_to_comp.end()) {
                        comp.v_to_e.push_back(std::vector<int>());
                        comp_to_graph.push_back(v2);
                        iter2 = graph_to_comp.insert(std::make_pair(v2,
                                    comp_to_graph.size()-1)).first;
                    }
                    comp.edges.push_back(Edge(iter1->second, iter2->second));
                    comp.v_to_e[iter1->second].push_back(comp.edges.size()-1);
                    comp.v_to_e[iter2->second].push_back(comp.edges.size()-1);
                    edge_stack.pop();
                    components.push_back(comp);
                    components_idx.push_back(comp_to_graph);
                }
            }
        }
    }
}

/* BFS for shortest path between start and end in subgraph. Returns Subgraph of
 * the shortest path with a vertex list for the path and, for each edge in the
 * path, sets path_no_clear.edges[edge] = true. Assumes that start and end are
 * connected in subgraph. */
void desres::msys::SSSR::get_subgraph_path(const GraphRepr& graph,
        const Subgraph& subgraph, int start, int end, Subgraph& path,
        Subgraph& path_no_clear) {
    /* Keep track of edge path in BFS, where the edge is represented as a 
     * (edge_idx, head_vertex_idx) pair */
    typedef std::vector<std::pair<int, int> > EdgeList;
    std::queue<EdgeList> BFS_queue;
    std::vector<bool> visited(graph.v_to_e.size(), false);
    std::pair<int, int> begin = std::make_pair(-1,start);
    BFS_queue.push(EdgeList(1,begin));
    visited[start] = true;
    while (BFS_queue.size() > 0) {
        const EdgeList& edges = BFS_queue.front();
        int v = edges[edges.size()-1].second;
        bool found = false;
        /* Loop through adjacent edges */
        for (unsigned i = 0; i < graph.v_to_e[v].size(); ++i) {
            int e = graph.v_to_e[v][i];
            if (!subgraph.edges[e]) continue;
            int other = graph.other(e,v);
            if (other == end) {
                BFS_queue.front().push_back(std::make_pair(e, other));
                found = true;
                break;
            }
            if (visited[other]) continue;
            visited[other] = true;
            BFS_queue.push(edges);
            BFS_queue.back().push_back(std::make_pair(e, other));
        }
        if (found) break;
        BFS_queue.pop();
    }
    if (BFS_queue.size() == 0)
        MSYS_FAIL("No path found between vertices in subgraph");
    /* Desired path is the list of edges in the first item of the queue */
    path.edges.clear();
    path.edges.resize(graph.edges.size(), false);
    path.vertex_list.resize(BFS_queue.front().size());
    path.vertex_list[0] = start;
    for (unsigned i = 1; i < BFS_queue.front().size(); ++i) {
        path.edges[BFS_queue.front()[i].first] = true;
        path.vertex_list[i] = BFS_queue.front()[i].second;
        if (!path_no_clear.edges[BFS_queue.front()[i].first])
            path_no_clear.edges[BFS_queue.front()[i].first] = true;
    }
}

/* BFS for shortest odd path(s) between start and end in graph with odd and even
 * edges. Returns Subgraph for a shortest odd path with the vertex list for the
 * path, or all shortest paths if return_multiple is true. Only returns shortest
 * path(s) with at most max_length edges, or all shortest path(s) if max_length
 * is -1. */
void desres::msys::SSSR::get_odd_path(const GraphRepr& graph,
        const Subgraph& odd_edges, int start, int end, bool return_multiple,
        int max_length, std::vector<Subgraph>& paths) {

    /* Represent a vertex as a (vertex-id, parity) pair, where true = odd */
    typedef std::pair<int, bool> Vertex;
    /* Keep track of edge path in BFS, where the edge is represented as a
     * (edge_idx, head_vertex) pair */
    typedef std::vector<std::pair<int, Vertex> > EdgeList;
    std::queue<EdgeList> BFS_queue;
    /* Keep track of the depth in the BFS of each visited vertex */
    std::vector<int> depths_odd(graph.v_to_e.size(), -1);
    std::vector<int> depths_even(graph.v_to_e.size(), -1);

    /* Start with an even vertex */
    std::pair<int, Vertex> begin(-1, Vertex(start, false));
    BFS_queue.push(EdgeList(1, begin));
    depths_even[start] = 0;
    int shortest_path_length = -1;
    while (BFS_queue.size() > 0) {
        const EdgeList& edges = BFS_queue.front();
        Vertex v = edges[edges.size()-1].second;
        /* Loop through adjacent edges */
        for (unsigned i = 0; i < graph.v_to_e[v.first].size(); ++i) {
            int e = graph.v_to_e[v.first][i];
            /* Switch vertex parity if e is an odd edge */
            bool odd = v.second ^ odd_edges.edges[e];
            int other = graph.other(e,v.first);
            /* The path is done only if we find end and parity is odd */
            if (other == end && odd) {
                Subgraph path(graph.edges.size(), false);
                path.vertex_list.resize(edges.size()+1);
                path.vertex_list[0] = start;
                for (unsigned i = 1; i < edges.size(); ++i) {
                    path.edges[edges[i].first] = true;
                    path.vertex_list[i] = edges[i].second.first;
                }
                path.edges[e] = true;
                path.vertex_list[path.vertex_list.size()-1] = end;
                paths.push_back(path);
                if (!return_multiple)
                    return;
                shortest_path_length = path.vertex_list.size() - 1;
            } else {
                int& depth = (odd ? depths_odd[other] : depths_even[other]);
                if (depth != -1 && depth < int(edges.size())) {
                    /* We have visited this vertex before at a lower depth */
                    continue;
                }
                depth = edges.size();
                if ((shortest_path_length == -1
                            || int(edges.size()) < shortest_path_length)
                        && (max_length == -1
                            || int(edges.size()) < max_length)) {
                    /* Only add vertex to BFS if path is less than max_length
                     * and length of shortest path */
                    BFS_queue.push(edges);
                    BFS_queue.back().push_back(std::make_pair(e,
                                Vertex(other, odd)));
                }
            }
        }
        BFS_queue.pop();
    }
}

/* Get a cycle basis for a graph. See F. Berger, P. Gritzmann, and S. de Vries,
 * "Minimum Cycle Bases for Network Graphs", Algorithmica 40, 2004, and
 * U. Bauer, "Minimum Cycle Basis Algorithms for the Chemistry Development
 * Toolkit", 2004.
 *
 * Assumes that graph is biconnected. Returns basis, a list of subgraphs each
 * representing a single ring; pivots, a corresponding list of edge indices 
 * where pivots[i] is an edge contained in basis[i] but not in basis[j] for any
 * j>i; and non_pivots, the list of non-pivot edge indices. Hence the 
 * basis-edge indicator matrix is upper-triangular if the columns are permuted
 * such that columns pivots[0],...,pivots[pivots.size()-1] come first. The 
 * returned basis is such that the first k basis rings are in fact minimal;
 * this value k is the unsigned return value. */
unsigned desres::msys::SSSR::get_cycle_basis(const GraphRepr& graph,
        std::deque<Subgraph>& basis, std::deque<int>& pivots,
        std::vector<int>& non_pivots) {

    basis.clear();
    pivots.clear();
    non_pivots.clear();
    Subgraph all(graph.edges.size(), true);
    Subgraph covered(graph.edges.size(), false);

    /* Step 1: Repeatedly pick an uncovered edge, add that edge to the front of
     * a pivot list, find the smallest cycle containing that edge, add that 
     * cycle to the front of the basis, and cover all edges in that cycle */
    Subgraph pivots_removed(graph.edges.size(), true);
    for (unsigned i = 0; i < graph.edges.size(); ++i) {
        if (covered.edges[i]) continue;
        pivots.push_front(i);
        pivots_removed.edges[i] = false;
        /* Find smallest cycle */
        basis.push_front(Subgraph());
        all.edges[i] = false;
        get_subgraph_path(graph, all, graph.edges[i].first,
                graph.edges[i].second, basis.front(), covered);
        basis.front().edges[i] = true;
        all.edges[i] = true;
    }

    /* Step 2: In the subgraph of non-pivot edges from step 1, form a spanning 
     * tree. For each non-pivot edge not in the tree, find the cycle consisting
     * of that edge and edges in the spanning tree, add that cycle to the back 
     * of the basis, and add that edge to the back of the pivot list */
    unsigned nfixed = pivots.size(); /* Save number of pivots from step 1 */
    Subgraph spanning_tree(graph.edges.size(), false);
    /* BFS for spanning tree */
    std::queue<int> BFS_queue;
    std::vector<bool> visited(graph.v_to_e.size(), false);
    BFS_queue.push(0);
    visited[0] = true;
    while (BFS_queue.size() > 0) {
        int v = BFS_queue.front();
        BFS_queue.pop();
        /* Loop through adjacent edges */
        for (unsigned i = 0; i < graph.v_to_e[v].size(); ++i) {
            int e = graph.v_to_e[v][i];
            if (!pivots_removed.edges[e]) continue;
            if (!visited[graph.other(e,v)]) {
                /* An unvisited edge is a tree edge, which is a non-pivot */
                spanning_tree.edges[e] = true;
                visited[graph.other(e,v)] = true;
                non_pivots.push_back(e);
                BFS_queue.push(graph.other(e,v));
            } else if (!spanning_tree.edges[e] && graph.other(e,v) < v) {
                /* Each non-tree edge is visited twice; include it once as a
                 * pivot */
                pivots.push_back(e);
            }
        }
    }
    /* For each new pivot edge, find cycle with that edge and edges in spanning
     * tree */
    for (unsigned i = nfixed; i < pivots.size(); ++i) {
        basis.push_back(Subgraph());
        get_subgraph_path(graph, spanning_tree, graph.edges[pivots[i]].first,
                graph.edges[pivots[i]].second, basis.back(), covered);
        basis.back().edges[pivots[i]] = true;
    }
    /* The cycles from step 1 are guaranteed minimal; the cycles from step 2
     * are not */
    return nfixed;
}

/* Minimize a cycle basis for a graph. See F. Berger, P. Gritzmann, and S. de
 * Vries, "Minimum Cycle Bases for Network Graphs", Algorithmica 40, 2004. 
 *
 * Given a ring basis of a graph, the first nfixed of which are minimal, and
 * lists of pivot and non-pivot edges such that the basis-edge indicator matrix
 * is upper-triangular when the columns are ordered pivots[0],...,
 * pivots[pivots.size()-1],..., returns a minimal basis with each non-minimal
 * ring replaced by a smaller ring. The input basis may be overwritten. */
void desres::msys::SSSR::minimize_cycle_basis(const GraphRepr& graph,
        std::deque<Subgraph>& basis, const std::deque<int>& pivots,
        unsigned nfixed, const std::vector<int>& non_pivots,
        std::vector<Subgraph>& min_basis) {

    min_basis.clear();
    for (unsigned i = 0; i < nfixed; ++i)
        min_basis.push_back(basis[i]);
    for (unsigned i = nfixed; i < pivots.size(); ++i) {
        /* Step 1: Construct a kernel vector U for the current basis with cycle
         * i removed */
        Subgraph U(graph.edges.size(), false);
        std::vector<bool> adj_vertices(graph.v_to_e.size(), false);
        /* Set U[pivots[i]] = 1 */
        U.edges[pivots[i]] = true;
        adj_vertices[graph.edges[pivots[i]].first] = true;
        adj_vertices[graph.edges[pivots[i]].second] = true;
        /* Solve upper-triangular system for U.edges[pivots[i-1]],...,
         * U.edges[pivots[0]] */
        for (int j = i-1; j >= 0; --j) {
            U.edges[pivots[j]] = basis[j].edges[pivots[i]];
            for (unsigned k = j+1; k < i; ++k)
                U.edges[pivots[j]] = U.edges[pivots[j]] ^
                    (U.edges[pivots[k]] & basis[j].edges[pivots[k]]);
            if (U.edges[pivots[j]]) {
                adj_vertices[graph.edges[pivots[j]].first] = true;
                adj_vertices[graph.edges[pivots[j]].second] = true;
            }
        }

        /* Step 2: Consider U as a set of odd edges. For all vertices v incident
         * to an edge in U, find the shortest odd cycle starting and ending at 
         * v and keep the shortest such cycle. */
        bool replaced = false;
        Subgraph shortest = basis[i];
        for (unsigned v = 0; v < graph.v_to_e.size(); ++v) {
            if (!adj_vertices[v]) continue;
            std::vector<Subgraph> current;
            get_odd_path(graph, U, v, v, false, -1, current);
            if (current.size() == 0)
                MSYS_FAIL("No odd path between vertices");
            current[0].vertex_list.erase(current[0].vertex_list.end()-1);
            if (current[0].vertex_list.size() < shortest.vertex_list.size()) {
                replaced = true;
                shortest = current[0];
            }
        }

        /* Step 3: If we have found a shorter cycle, replace the basis cycle 
         * with this one and do Gaussian elimination */
        if (replaced) {
            min_basis.push_back(shortest);
            for (unsigned j = 0; j < pivots.size(); ++j) {
                if (j < i && shortest.edges[pivots[j]]) {
                    for (unsigned k = j; k < pivots.size(); ++k)
                        shortest.edges[pivots[k]] = shortest.edges[pivots[k]]
                            ^ basis[j].edges[pivots[k]];
                    for (unsigned k = 0; k < non_pivots.size(); ++k)
                        shortest.edges[non_pivots[k]]
                            = shortest.edges[non_pivots[k]]
                            ^ basis[j].edges[non_pivots[k]];
                }
                basis[i].edges[pivots[j]] = shortest.edges[pivots[j]];
            }
            for (unsigned j = 0; j < non_pivots.size(); ++j)
                basis[i].edges[non_pivots[j]] = shortest.edges[non_pivots[j]];
        } else
            min_basis.push_back(basis[i]);
    }
}

/* Get all relevant cycles from a minimal cycle basis. See U. Bauer, "Minimum
 * Cycle Basis Algorithms for the Chemistry Development Toolkit", 2004. The
 * input min_basis may be overwritten. */
void desres::msys::SSSR::get_relevant_cycles(const GraphRepr& graph,
        std::vector<Subgraph>& min_basis, const std::deque<int>& pivots,
        const std::vector<int>& non_pivots,
        std::set<Subgraph>& relevant_cycles) {

    /* Step 1: Given the cycle-edge indicator matrix, compute a kernel vector
     * for this matrix with each cycle removed. */
    std::vector<Subgraph> kernels(pivots.size(), Subgraph(graph.edges.size(),
                false));
    /* Initialize kernel matrix (with the kernel vectors as columns) to the
     * identity on the pivot rows and 0 on the non-pivot rows */
    for (unsigned i = 0; i < pivots.size(); ++i)
        kernels[i].edges[pivots[i]] = true;
    /* Use Gaussian elimination to invert pivot-columns-submatrix of cycle-edge
     * indicator matrix, simultaneously applying operations to kernel matrix to
     * obtain inverse */
    for (unsigned i = 0; i < pivots.size(); ++i) {
        /* Swap row i with some row j > i if necessary */
        if (!min_basis[i].edges[pivots[i]]) {
            unsigned j;
            for (j = i+1; j < pivots.size(); ++j) {
                if (min_basis[j].edges[pivots[i]])
                    break;
            }
            if (j == pivots.size())
                MSYS_FAIL("Min basis indicator matrix is singular");
            std::swap(min_basis[i].edges, min_basis[j].edges);
            for (unsigned k = 0; k < pivots.size(); ++k) {
                bool tmp = kernels[k].edges[pivots[i]];
                kernels[k].edges[pivots[i]] = kernels[k].edges[pivots[j]];
                kernels[k].edges[pivots[j]] = tmp;
            }
        }
        /* Eliminate 1's in column i, rows j > i */
        for (unsigned j = i+1; j < pivots.size(); ++j) {
            if (min_basis[j].edges[pivots[i]]) {
                for (unsigned k = 0; k < pivots.size(); ++k) {
                    if (k >= i)
                        min_basis[j].edges[pivots[k]]
                            = min_basis[j].edges[pivots[k]]
                            ^ min_basis[i].edges[pivots[k]];
                    kernels[k].edges[pivots[j]] = kernels[k].edges[pivots[j]]
                        ^ kernels[k].edges[pivots[i]];
                }
            }
        }
    }
    for (unsigned i = pivots.size()-1; i > 0; --i) {
        /* Eliminate 1's in column i, rows j < i */
        for (unsigned j = 0; j < i; ++j) {
            if (min_basis[j].edges[pivots[i]]) {
                for (unsigned k = 0; k < pivots.size(); ++k) {
                    if (k >= i)
                        min_basis[j].edges[pivots[k]]
                            = min_basis[j].edges[pivots[k]]
                            ^ min_basis[i].edges[pivots[k]];
                    kernels[k].edges[pivots[j]] = kernels[k].edges[pivots[j]]
                        ^ kernels[k].edges[pivots[i]];
                }
            }
        }
    }

    /* Step 2: Consider each kernel vector U as a subgraph of odd edges. For
     * each vertex v incident to an edge in U, find all (shortest) odd cycles
     * starting and ending at v with length equal to that of the original basis
     * vector corresponding to U. */
    std::vector<Subgraph> paths;
    for (unsigned i = 0; i < pivots.size(); ++i) {
        std::set<int> checked;
        for (unsigned j = 0; j < kernels[i].edges.size(); ++j) {
            if (!kernels[i].edges[j]) continue;
            int v1 = graph.edges[j].first;
            int v2 = graph.edges[j].second;
            if (checked.find(v1) == checked.end()) {
                checked.insert(v1);
                get_odd_path(graph, kernels[i], v1, v1, true,
                        min_basis[i].vertex_list.size(), paths);
            }
            if (checked.find(v2) == checked.end()) {
                checked.insert(v2);
                get_odd_path(graph, kernels[i], v2, v2, true,
                        min_basis[i].vertex_list.size(), paths);
            }
        }
    }

    /* Step 3: Filter the list of cycles for duplicates */
    relevant_cycles.clear();
    for (unsigned i = 0; i < paths.size(); ++i) {
        /* Reorder each cycle such that the first vertex has lowest index and
         * the smaller-indexed neighbor of this vertex is second */
        std::vector<int> path = paths[i].vertex_list;
        path.erase(path.end()-1);
        int min = graph.v_to_e.size();
        int min_ind = -1;
        for (unsigned j = 0; j < path.size(); ++j) {
            if (path[j] < min) {
                min = path[j];
                min_ind = j;
            }
        }
        int right = (min_ind == int(path.size()) ? 0 : min_ind + 1);
        int left = (min_ind == 0 ? path.size()-1 : min_ind - 1);
        paths[i].vertex_list.clear();
        paths[i].vertex_list.reserve(path.size());
        paths[i].vertex_list.push_back(min);
        if (path[right] < path[left]) {
            for (int j = right; j != min_ind; j = (j+1) % path.size())
                paths[i].vertex_list.push_back(path[j]);
        } else {
            for (int j = left; j != min_ind; j = (j == 0 ? path.size()-1 : j-1))
                paths[i].vertex_list.push_back(path[j]);
        }

        relevant_cycles.insert(paths[i]);
    }
}

desres::msys::MultiIdList 
desres::msys::GetSSSR(SystemPtr mol, IdList const& atoms, 
        bool all_relevant) {

    MultiIdList sssr;

    /* keep atoms that are not virtuals or transistion metals */
    bondedVirtualsAndMetalsFilter keep(mol);

    /* Create GraphRepr with reindexed atoms (vertices) and bonds (edges) */
    std::vector<int> atom_idx_map(mol->maxAtomId(), -1);
    for (unsigned i = 0; i < atoms.size(); ++i) {
        if (mol->atom(atoms[i]).atomic_number >= 1)
            atom_idx_map[atoms[i]] = i;
    }
    GraphRepr graph;
    graph.v_to_e.resize(atoms.size(), std::vector<int>());
    for (unsigned i = 0; i < atoms.size(); ++i) {
        /* If i is a pseudo atom or metal, graph.v_to_e[i] remains empty. These graph
         * vertex indices are ignored by get_biconnected_components. */
        if (!keep(mol->atom(atoms[i]))) continue;
        IdList bonded = mol->filteredBondedAtoms(atoms[i],keep);
        for (unsigned j = 0; j < bonded.size(); ++j) {
            if (atoms[i] < bonded[j] && atom_idx_map[bonded[j]] != -1) {
                graph.edges.push_back(Edge(i, atom_idx_map[bonded[j]]));
                graph.v_to_e[i].push_back(graph.edges.size()-1);
                graph.v_to_e[atom_idx_map[bonded[j]]].push_back(
                        graph.edges.size()-1);
            }
        }
    }

    /* Get biconnected components */
    std::vector<GraphRepr> components;
    std::vector<std::vector<int> > components_idx;
    get_biconnected_components(graph, components, components_idx);

    /* Get SSSR for each biconnected component */
    for (unsigned i = 0; i < components.size(); ++i) {
        if (components[i].v_to_e.size() <= 2) continue;
        std::deque<Subgraph> basis;
        std::deque<int> pivots;
        std::vector<int> non_pivots;
        unsigned nfixed = get_cycle_basis(components[i], basis, pivots, 
                non_pivots);
        std::vector<Subgraph> min_basis;
        minimize_cycle_basis(components[i], basis, pivots, nfixed, non_pivots,
                min_basis);
        if (!all_relevant) {
            for (unsigned j = 0; j < min_basis.size(); ++j) {
                sssr.push_back(IdList());
                for (unsigned k = 0; k < min_basis[j].vertex_list.size(); ++k) {
                    /* Map back to original atom IDs */
                    int v = min_basis[j].vertex_list[k];
                    sssr[sssr.size()-1].push_back(atoms[components_idx[i][v]]);
                }
            }
        } else {
            std::set<Subgraph> relevant_rings;
            get_relevant_cycles(components[i], min_basis, pivots, non_pivots,
                    relevant_rings);
            for (std::set<Subgraph>::iterator iter = relevant_rings.begin();
                    iter != relevant_rings.end(); ++iter) {
                sssr.push_back(IdList());
                for (unsigned k = 0; k < iter->vertex_list.size(); ++k) {
                    /* Map back to original atom IDs */
                    int v = iter->vertex_list[k];
                    sssr[sssr.size()-1].push_back(atoms[components_idx[i][v]]);
                }
            }
        }
    }
    return sssr;
}

desres::msys::MultiIdList desres::msys::FusedRingSystems(SystemPtr mol,
        MultiIdList const& rings) {

    MultiIdList bond_to_rings(mol->maxBondId());
    for (unsigned i = 0; i < rings.size(); ++i) {
        for (unsigned j = 0; j < rings[i].size(); ++j) {
            Id bond = mol->findBond(rings[i][j],
                    rings[i][(j+1)%rings[i].size()]);
            if (bond == msys::BadId)
                MSYS_FAIL("Ring bond not found in system");
            bond_to_rings[bond].push_back(i);
        }
    }

    std::vector<bool> processed_bonds(mol->maxBondId(), false);
    MultiIdList ring_systems;
    BOOST_FOREACH(Id bond, mol->bonds()) {
        if (processed_bonds[bond]) continue;
        processed_bonds[bond] = true;
        if (bond_to_rings[bond].size() == 0) continue;
        /* Get the ring system containing this bond */
        std::set<Id> bond_set;
        std::set<Id> ring_set;
        std::queue<Id> unprocessed_bonds;
        bond_set.insert(bond);
        unprocessed_bonds.push(bond);
        while (unprocessed_bonds.size() > 0) {
            Id front = unprocessed_bonds.front();
            unprocessed_bonds.pop();
            /* Loop through all potentially aromatic SSSR rings containing
             * bond 'front' */
            BOOST_FOREACH(Id ring, bond_to_rings[front]) {
                ring_set.insert(ring);
                for (unsigned i = 0; i < rings[ring].size(); ++i) {
                    Id ring_bond = mol->findBond(rings[ring][i],
                            rings[ring][(i+1)%rings[ring].size()]);
                    /* If ring bond is new, add to unprocessed bond queue */
                    if (bond_set.insert(ring_bond).second) {
                        unprocessed_bonds.push(ring_bond);
                        processed_bonds[ring_bond] = true;
                    }
                }
            }
        }
        /* Add this ring system */
        ring_systems.push_back(IdList(ring_set.begin(), ring_set.end()));
    }
    return ring_systems;
}
