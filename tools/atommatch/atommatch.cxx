#include "atommatch.hxx"
#include "atommatch_helpers.hxx"
#include <queue>
#include <stack>
#include <limits>
#include <cmath>
#include <Eigen/Dense>

using namespace desres::msys;
using namespace desres::msys::atommatch;
using namespace desres::msys::atommatch::helpers;

/* Create a graph from a set of connected atoms in a molecule.
 *
 * g is the created graph object and atom_idx is the atom ID mapping such
 * that node i in g maps to atom atom_idx[i] in mol.
 */
void helpers::init_graph(SystemPtr mol, std::vector<int>& atom_idx,
        GraphRepr& g) {
    atom_idx.clear();
    atom_idx.resize(mol->maxAtomId(), -1);
    IdList atoms = mol->atoms();
    for (unsigned i = 0; i < atoms.size(); ++i)
        atom_idx[atoms[i]] = i;
    g.v_to_e.clear();
    g.edges.clear();
    g.v_to_e.resize(atoms.size(), std::vector<int>());
    for (unsigned i = 0; i < atoms.size(); ++i) {
        IdList bonded = mol->bondedAtoms(atoms[i]);
        BOOST_FOREACH(Id other, bonded) {
            if (atoms[i] < other) {
                g.edges.push_back(Edge(i, atom_idx[other]));
                g.v_to_e[i].push_back(g.edges.size() - 1);
                g.v_to_e[atom_idx[other]].push_back(g.edges.size() - 1);
            }
        }
    }
}

/* Partition a graph g into biconnected components.
 *
 * components[i] is the graph of the i'th biconnected component of g, where the
 * components are ordered in decreasing size. components_idx[i] is the ID
 * mapping such that node j of components[i] corresponds to node
 * components_idx[i][j] of g.
 *
 * tree is a graph of biconnected components, where node i of tree represents
 * the i'th biconnected component of g and nodes i and j of tree are connected
 * by an edge iff components i and j of g share a node. Note that tree is not
 * quite a tree but a 'tree-like' graph, in the following sense: If a single
 * node of tree is picked as the root and all other nodes are assigned a depth
 * value representing the edge distance to the root, then each node at depth
 * d > 0 is connected to a single parent node at depth d-1, and the remaining
 * neighbors of the node are at depth d or d+1.
 *
 * tree_edges[i][j] is a list of pairs representing components (other than
 * component i) of g that share node j of components[i]; the first index in each
 * pair is the ID of the other component and the second is the ID of the shared
 * node in the other component.
 */
void helpers::partition_graph_biconnected(const GraphRepr& g,
        std::vector<GraphRepr>& components,
        std::vector<std::vector<int> >& components_idx, GraphRepr& tree,
        TreeEdges& tree_edges) {
    /* Use msys routine to get biconnected components */
    std::vector<GraphRepr> tmp_components;
    std::vector<std::vector<int> > tmp_components_idx;
    SSSR::get_biconnected_components(g, tmp_components,
            tmp_components_idx);
    std::map<int, std::vector<int> > by_size;
    for (unsigned i = 0; i < tmp_components_idx.size(); ++i)
        by_size[tmp_components_idx[i].size()].push_back(i);
    /* Sort components by size; this is only important in that if RMSD
     * alignment is used to select a match from multiple best matches, then the
     * heuristic algorithm to do so will try to align the largest matching
     * biconnected components first. */
    components.clear();
    components_idx.clear();
    for (std::map<int, std::vector<int> >::reverse_iterator
            iter = by_size.rbegin(); iter != by_size.rend(); ++iter) {
        BOOST_FOREACH(int i, iter->second) {
            components.push_back(tmp_components[i]);
            components_idx.push_back(tmp_components_idx[i]);
        }
    }
    /* Construct tree of biconnected components */
    tree.v_to_e.clear();
    tree.v_to_e.resize(components.size());
    tree.edges.clear();
    tree_edges.clear();
    tree_edges.resize(components.size());
    for (unsigned i = 0; i < tree_edges.size(); ++i)
        tree_edges[i].resize(components_idx[i].size());
    std::vector<std::vector<std::pair<int, int> > > idx(g.v_to_e.size());
    for (unsigned i = 0; i < components_idx.size(); ++i) {
        for (unsigned j = 0; j < components_idx[i].size(); ++j) {
            const std::vector<std::pair<int, int> >& prev
                = idx[components_idx[i][j]];
            for (unsigned k = 0; k < prev.size(); ++k) {
                tree.edges.push_back(Edge(i, prev[k].first));
                tree.v_to_e[i].push_back(tree.edges.size()-1);
                tree.v_to_e[prev[k].first].push_back(tree.edges.size()-1);
                tree_edges[i][j].push_back(prev[k]);
                tree_edges[prev[k].first][prev[k].second].push_back(
                        std::make_pair(i,j));
            }
            idx[components_idx[i][j]].push_back(std::make_pair(i,j));
        }
    }
}

/* Partition a graph into bi-edge-connected components.
 *
 * This function is currently not used; keep it around in case we switch the
 * method later.
 */
#if 0
void helpers::partition_graph_biedgeconnected(const GraphRepr& g,
        std::vector<GraphRepr>& components,
        std::vector<std::vector<int> >& components_idx, GraphRepr& tree,
        std::vector<std::vector<std::pair<int, int> > >& tree_edges) {
    /* Tarjan algorithm for all bridge edges
     * See Tarjan, "A note on finding the bridges of a graph", Information
     * Processing Letters 2(1974) 160-161. */
    /* Label nodes in pre-order (rather than post-order as in Tarjan paper) */
    std::vector<bool> in_tree(g.edges.size(), false);
    std::vector<int> v_to_label(g.v_to_e.size(), -1);
    std::vector<int> label_to_v;
    label_to_v.reserve(g.v_to_e.size());
    std::stack<std::pair<int, int> > S;
    S.push(std::make_pair(0,0));
    label_to_v.push_back(0);
    v_to_label[0] = 0;
    while (true) {
        while (S.top().second >= int(g.v_to_e[S.top().first].size())) {
            S.pop();
            if (S.size() == 0) break;
            S.top().second += 1;
        }
        if (S.size() == 0) break;
        int edge = g.v_to_e[S.top().first][S.top().second];
        int other = g.other(edge, S.top().first);
        if (v_to_label[other] != -1)
            S.top().second += 1;
        else {
            S.push(std::make_pair(other,0));
            v_to_label[other] = label_to_v.size();
            label_to_v.push_back(other);
            in_tree[edge] = true;
        }
    }
    if (label_to_v.size() != g.v_to_e.size())
        MSYS_FAIL("Input molecule is not connected");
    /* Compute ND, L, and H as in Tarjan paper */
    std::vector<int> ND(g.v_to_e.size(), 0);
    std::vector<int> low(g.v_to_e.size(), -1);
    std::vector<int> high(g.v_to_e.size(), -1);
    std::vector<bool> is_bridge(g.edges.size(), false);
    for (int i = label_to_v.size() - 1; i >= 0; --i) {
        int nd_sum = 1;
        int low_min = i;
        int high_max = i;
        int back_edge = -1;
        BOOST_FOREACH(int edge, g.v_to_e[label_to_v[i]]) {
            int other = g.other(edge, label_to_v[i]);
            if (in_tree[edge] && v_to_label[other] > i) {
                /* Forward tree edge */
                nd_sum += ND[other];
                low_min = (low_min < low[other] ? low_min : low[other]);
                high_max = (high_max > high[other] ? high_max : high[other]);
            } else if (in_tree[edge]) {
                /* Backward tree edge */
                back_edge = edge;
            } else {
                /* Non-tree edge */
                low_min = (low_min < v_to_label[other]
                        ? low_min : v_to_label[other]);
                high_max = (high_max > v_to_label[other]
                        ? high_max : v_to_label[other]);
            }
        }
        ND[label_to_v[i]] = nd_sum;
        low[label_to_v[i]] = low_min;
        high[label_to_v[i]] = high_max;
        if (low_min == i && high_max < i + nd_sum && back_edge != -1)
            /* Tree edge pointing to this vertex is a bridge */
            is_bridge[back_edge] = true;
    }

    /* Get graphs of 2-edge-connected components and tree of such components */
    components.clear();
    components_idx.clear();
    tree.v_to_e.clear();
    tree.edges.clear();
    tree_edges.clear();
    tree_edges.resize(g.v_to_e.size());
    std::vector<int> component_id(g.v_to_e.size(), -1);
    std::vector<int> lookup(g.v_to_e.size(), -1);
    for (unsigned i = 0; i < g.v_to_e.size(); ++i) {
        if (component_id[i] != -1) continue;
        /* Create new component graph and tree node */
        components.push_back(GraphRepr());
        components_idx.push_back(std::vector<int>());
        tree.v_to_e.push_back(std::vector<int>());
        /* BFS for this component */
        std::queue<int> Q;
        Q.push(i);
        component_id[i] = components.size()-1;
        components.back().v_to_e.push_back(std::vector<int>());
        components_idx.back().push_back(i);
        lookup[i] = 0;
        while (Q.size() > 0) {
            int top = Q.front();
            Q.pop();
            BOOST_FOREACH(int edge, g.v_to_e[top]) {
                int other = g.other(edge, top);
                if (is_bridge[edge] && component_id[other] != -1) {
                    /* Edge to already constructed component--add tree edge */
                    tree.edges.push_back(Edge(tree.v_to_e.size()-1,
                                component_id[other]));
                    tree.v_to_e.back().push_back(tree.edges.size()-1);
                    tree.v_to_e[component_id[other]].push_back(
                            tree.edges.size()-1);
                    tree_edges[top].push_back(std::make_pair(
                                component_id[other], lookup[other]));
                    tree_edges[other].push_back(std::make_pair(
                                tree.v_to_e.size()-1, lookup[top]));
                } else if (!is_bridge[edge] && component_id[other] != -1) {
                    /* Edge to already visited node in this component */
                    if (lookup[top] < lookup[other]) {
                        components.back().edges.push_back(Edge(lookup[top],
                                    lookup[other]));
                        components.back().v_to_e[lookup[top]].push_back(components.back().edges.size()-1);
                        components.back().v_to_e[lookup[other]].push_back(components.back().edges.size()-1);
                    }
                } else if (!is_bridge[edge]) {
                    /* Edge to new node in this component--add component node
                     * and edge, and push node on BFS queue */
                    component_id[other] = components.size()-1;
                    Q.push(other);
                    components.back().v_to_e.push_back(std::vector<int>());
                    components.back().edges.push_back(Edge(lookup[top],
                                components.back().v_to_e.size()-1));
                    components.back().v_to_e.back().push_back(
                            components.back().edges.size()-1);
                    components.back().v_to_e[lookup[top]].push_back(
                            components.back().edges.size()-1);
                    components_idx.back().push_back(other);
                    lookup[other] = components.back().v_to_e.size()-1;
                }
            }
        }
    }
}
#endif

/* Helper routine for isomorphisms function.
 *
 * Determines whether node v1 can be matched to node v2 given the matches
 * found so far.
 */
bool match_node(const GraphRepr& g1, int v1, const GraphRepr& g2, int v2,
        const std::vector<int>& G1toG2, const std::vector<int>& G2toG1) {
    /* Check node degrees are the same */
    if (g1.v_to_e[v1].size() != g2.v_to_e[v2].size())
        return false;
    BOOST_FOREACH(int edge1, g1.v_to_e[v1]) {
        int hnbr = G1toG2[g1.other(edge1, v1)];
        if (hnbr == -1)
            continue;
        /* Check that matched neighbors of g map to neighbors of h */
        bool is_nbr = false;
        BOOST_FOREACH(int edge2, g2.v_to_e[v2]) {
            if (g2.other(edge2, v2) == hnbr)
                is_nbr = true;
        }
        if (!is_nbr)
            return false;
    }
    BOOST_FOREACH(int edge2, g2.v_to_e[v2]) {
        int gnbr = G2toG1[g2.other(edge2, v2)];
        if (gnbr == -1)
            continue;
        /* Check that matched neighbors of h map to neighbors of g */
        bool is_nbr = false;
        BOOST_FOREACH(int edge1, g1.v_to_e[v1]) {
            if (g1.other(edge1, v1) == gnbr)
                is_nbr = true;
        }
        if (!is_nbr)
            return false;
   }
   return true;
}

/* Find and score all isomorphisms between two graphs g1 and g2.
 *
 * Isomorphisms are determined only by bond structure; they are not affected by
 * atomic number or other atom/bond properties. Score for each isomorphism is
 * given by score_fct input.
 *
 * comp_idx1, atom_idx1, mol1, comp_idx2, atom_idx2, and mol2 inputs are only
 * used to map node IDs in g1 and g2 back to original atoms, for use in
 * score_fct.
 *
 * isos[i] is the i'th isomorphism, where each node j of g1 maps to
 * node isos[i].perm[j] of g2, and isos[i].score is the score of the
 * isomorphism.
 */
void helpers::isomorphisms(const GraphRepr& g1,
        const std::vector<int>& comp_idx1,
        const std::vector<int>& atom_idx1, SystemPtr mol1,
        const GraphRepr& g2, const std::vector<int>& comp_idx2,
        const std::vector<int>& atom_idx2, SystemPtr mol2,
        ScoreFctPtr score_fct, Isomorphisms& isos) {
    isos.clear();
    if (g1.edges.size() != g2.edges.size()
            || g1.v_to_e.size() != g2.v_to_e.size())
        return;
    int size = g1.v_to_e.size();
    
    /* Nodes of g1 in BFS order from node 0 */
    std::vector<bool> visited(size, false);
    std::queue<int> Q;
    std::vector<int> nodes_bfs;
    nodes_bfs.reserve(size);
    Q.push(0);
    nodes_bfs.push_back(0);
    visited[0] = true;
    while (Q.size() > 0) {
        int front = Q.front();
        Q.pop();
        BOOST_FOREACH(int edge, g1.v_to_e[front]) {
            if (!visited[g1.other(edge, front)]) {
                Q.push(g1.other(edge, front));
                nodes_bfs.push_back(g1.other(edge, front));
                visited[g1.other(edge, front)] = true;
            }
        }
    }
    if (int(nodes_bfs.size()) != size)
        MSYS_FAIL("Graph isomorphism error---check that graph is connected");

    /* Maintain the following: For all but top element of matches, matches[i] of
     * H matches nodes_bfs[i] of G. Top element of matches indicates the node in
     * H we are currently trying to match to nodes_bfs[i] of G. */
    std::stack<int> matches;
    matches.push(0);
    std::vector<int> G1toG2(size, -1);
    std::vector<int> G2toG1(size, -1);
    while (matches.size() > 0) {
        int v1 = nodes_bfs[matches.size() - 1];
        int v2 = matches.top();
        if (G2toG1[v2] == -1 && match_node(g1, v1, g2, v2, G1toG2, G2toG1)) {
            /* Matches; try to match next node of G in nodes_bfs */
            G1toG2[v1] = v2;
            G2toG1[v2] = v1;
            matches.push(0);
        } else {
            /* Does not match; increment top element of matches */
            matches.pop();
            matches.push(v2+1);
        }
        if (int(matches.size()) == size + 1) {
            /* Have matched all nodes in g1 */
            isos.push_back(Isomorphism());
            isos.back().perm = G1toG2;
            IdList atoms1(size);
            IdList atoms2(size);
            for (int i = 0; i < size; ++i) {
                atoms1[i] = atom_idx1[comp_idx1[i]];
                atoms2[i] = atom_idx2[comp_idx2[G1toG2[i]]];
            }
            isos.back().score = score_fct->apply(mol1, atoms1, mol2, atoms2);

            /* Top element is the placeholder 0; remove it */
            matches.pop();

            /* Undo top match and increment */
            v2 = matches.top();
            G1toG2[G2toG1[v2]] = -1;
            G2toG1[v2] = -1;
            matches.pop();
            matches.push(v2+1);
        }
        while (matches.top() == size) {
            /* Top element is size of g2; remove it to return to previous node
             * of g1 in nodes_bfs */
            matches.pop();
            if (matches.size() == 0) break;
            /* Undo top match and increment */
            v2 = matches.top();
            G1toG2[G2toG1[v2]] = -1;
            G2toG1[v2] = -1;
            matches.pop();
            matches.push(v2+1);
        }
    }
}

/* Given a graph tree and a node root, determine the edge distance of each node
 * from root.
 *
 * levels[i] is the list of all nodes at distance i from root. level_map[i] is
 * the distance of node i from root (so levels[level_map[i]] contains i). If
 * screen is true, then all nodes in tree with node index less than root, as
 * well as their descendants in tree, are ignored and assigned -1 in level_map.
 */
void helpers::get_levels(const GraphRepr& tree, int root,
        std::vector<std::vector<int> >& levels, std::vector<int>& level_map,
        bool screen) {
    levels.clear();
    level_map.clear();
    level_map.resize(tree.v_to_e.size(), -1);
    level_map[root] = 0;
    levels.push_back(std::vector<int>(1, root));
    while(true) {
        std::vector<int> next;
        BOOST_FOREACH(int node, levels[levels.size()-1]) {
            BOOST_FOREACH(int edge, tree.v_to_e[node]) {
                int other = tree.other(edge, node);
                if (screen && other < root)
                    continue;
                if (level_map[other] == -1) {
                    next.push_back(other);
                    level_map[other] = levels.size();
                }
            }
        }
        if (next.size() > 0)
            levels.push_back(next);
        else
            break;
    }
}

/* Get all matchings of [0,1,...,size1-1] with [0,1,...,size2-1], where
 * elements of either list can be unmatched.
 *
 * perms[i] is the i'th matching, where element j in the first list is matched
 * to element perms[i][j] in the second list if perms[i][j] >= 0, or element j
 * in the first list is unmatched if perms[i][j] == -1.
 *
 * E.g., if size1 == 3 and size2 == 2, then this returns the following matches:
 * (-1,-1,-1)
 * (-1,-1,0)
 * (-1,-1,1)
 * (-1,0,-1)
 * (-1,1,-1)
 * (0,-1,-1)
 * (1,-1,-1)
 * (0,1,-1)
 * (1,0,-1)
 * (0,-1,1)
 * (1,-1,0)
 * (-1,0,1)
 * (-1,1,0)
 */
void helpers::get_perms(int size1, int size2,
        std::vector<std::vector<int> >& perms) {
    perms.clear();
    if (size1 == 0) {
        perms.push_back(std::vector<int>());
        return;
    }
    static std::map<std::pair<int, int>, std::vector<std::vector<int> > > Cache;
    std::map<std::pair<int, int>, std::vector<std::vector<int> > >::iterator
        iter = Cache.find(std::make_pair(size1, size2));
    if (iter != Cache.end()) {
        perms = iter->second;
        return;
    }
    std::vector<int> p(size1, -1);
    std::map<int, int> counts;
    counts[-1] = size1;
    while (true) {
        /* Check that no two elements of set A map to the same element of
         * set B */
        bool valid = true;
        for (std::map<int, int>::iterator iter = counts.begin();
                iter != counts.end(); ++iter) {
            if (iter->first != -1 && iter->second > 1)
                valid = false;
        }
        if (valid)
            perms.push_back(p);
        /* Increment to next mapping */
        int ind = 0;
        if (counts[p[ind]] == 1)
            counts.erase(p[ind]);
        else
            --counts[p[ind]];
        ++p[ind];
        while (p[ind] == size2) {
            p[ind] = -1;
            ++counts[-1];
            ++ind;
            if (ind == size1) break;
            if (counts[p[ind]] == 1)
                counts.erase(p[ind]);
            else
                --counts[p[ind]];
            ++p[ind];
        }
        if (ind == size1) break;
        ++counts[p[ind]];
    }
    Cache.insert(std::make_pair(std::make_pair(size1, size2), perms));
}

/* Helper function used in dp_update.
 *
 * Given a list of options of length n (where each option represents a list of
 * isomorphisms of biconnected components) and m possible isomorphisms of
 * biconnected component v1 of mol1 with v2 of mol2, updates options into a
 * list of m*n new options by appending each of the m isomorphisms of v1 and v2
 * to each original option. */
void helpers::options_append(std::vector<DPNodes>& options, int v1, int v2,
        const std::vector<int>& option_inds) {
    std::vector<DPNodes> return_options;
    BOOST_FOREACH(int ind, option_inds) {
        BOOST_FOREACH(const DPNodes& option, options) {
            return_options.push_back(option);
            return_options.back().push_back(DPNode(v1, v2, ind));
        }
    }
    options = return_options;
}

/* Core step of dynamic programming algorithm. (See max_match for details.)
 *
 * Update best_matches[v1][v2], assuming best_matches for descendant components
 * of v1 and v2 in the biconnected component trees of mol1 and mol2 have
 * already been updated.
 *
 * If keep_options is false, then this routine only keeps track of the maximum
 * isomorphism scores and sizes, without keeping track of what the isomorphisms
 * actually are that can attain these scores/sizes.
 */
void helpers::dp_update(int v1, int v2, DPMat& best_matches,
        const std::vector<std::vector<std::pair<int, int> > >& tree_edges1, 
        const std::vector<int>& level_map1,
        const std::vector<std::vector<std::pair<int, int> > >& tree_edges2,
        const std::vector<int>& level_map2,
        const IsoMat& iso_mat, bool keep_options) {
    best_matches[v1][v2].clear();
    best_matches[v1][v2].resize(iso_mat[v1][v2].size());
    for (unsigned iso_ind = 0; iso_ind < iso_mat[v1][v2].size(); ++iso_ind) {
        /* Update best_matches[v1][v2][iso_ind] for the iso_ind'th isomorphism
         * of components v1 and v2 */
        const std::vector<int>& perm = iso_mat[v1][v2][iso_ind].perm;
        std::vector<DPNodes>& options = best_matches[v1][v2][iso_ind].options;
        options.push_back(DPNodes());
        double& final_score = best_matches[v1][v2][iso_ind].score;
        final_score = iso_mat[v1][v2][iso_ind].score;
        int& final_size = best_matches[v1][v2][iso_ind].size;
        final_size = iso_mat[v1][v2][iso_ind].perm.size();
        /* For each matching pair (i,perm[i]) of atoms in this isomorphism, find
         * best matches between child components of v1 sharing atom i and child
         * components of v2 sharing atom perm[i] */
        for (unsigned i = 0; i < perm.size(); ++i) {
            std::vector<std::pair<int, int> > te1 = tree_edges1[i];
            std::vector<std::pair<int, int> > te2 = tree_edges2[perm[i]];
            for (std::vector<std::pair<int, int> >::iterator iter
                    = te1.begin(); iter != te1.end(); ) {
                if (level_map1[iter->first] <= level_map1[v1])
                    iter = te1.erase(iter);
                else
                    ++iter;
            }
            for (std::vector<std::pair<int, int> >::iterator iter
                    = te2.begin(); iter != te2.end(); ) {
                if (level_map2[iter->first] <= level_map2[v2])
                    iter = te2.erase(iter);
                else
                    ++iter;
            }
            /* Consider all possible matchings of child components of v1 sharing
             * atom i with child components of v2 sharing atom perm[i] */
            std::vector<std::vector<int> > comp_perms;
            get_perms(te1.size(), te2.size(), comp_perms);
            double max_score = 0;
            int max_size = 0;
            std::vector<const std::vector<int>* > max_perms;
            std::vector<std::vector<std::vector<int> > > max_option_inds;
            BOOST_FOREACH(const std::vector<int>& comp_perm, comp_perms) {
                double score = 0;
                int size = 0;
                std::vector<std::vector<int> > option_inds;
                bool valid = true;
                for (unsigned k = 0; k < comp_perm.size(); ++k) {
                    /* Determine best isomorphism of child component k of v1
                     * with child component comp_perm[k] of v2 that maps atom
                     * i of v1 to atom perm[i] of v2 */
                    option_inds.push_back(std::vector<int>());
                    if (comp_perm[k] == -1) continue;
                    double tmp_max_score = 0;
                    int tmp_max_size = 0;
                    /* All isomorphisms of child component k of v1 with child
                     * component comp_perm[k] of v2 */
                    const std::vector<DPBest>& isos
                        = best_matches[te1[k].first][te2[comp_perm[k]].first];
                    for (unsigned j = 0; j < isos.size(); ++j) {
                        const Isomorphism& iso
                            = iso_mat[te1[k].first][te2[comp_perm[k]].first][j];
                        /* Check that this isomorphism maps atom i of v1 to atom
                         * perm[i] of v2 */
                        if (iso.perm[te1[k].second]
                                == te2[comp_perm[k]].second) {
                            if (isos[j].score > tmp_max_score ||
                                    (isos[j].score == tmp_max_score
                                     && isos[j].size > tmp_max_size)) {
                                tmp_max_score = isos[j].score;
                                tmp_max_size = isos[j].size;
                                option_inds.back() = std::vector<int>(1,j);
                            } else if (isos[j].score == tmp_max_score
                                    && isos[j].size == tmp_max_size)
                                option_inds.back().push_back(j);
                        }
                    }
                    /* option_inds now contains all isomorphisms of child
                     * component k of v1 with child component comp_perm[k] of
                     * v2 that gives the best subtree score/size */
                    if (option_inds.back().size() == 0) {
                        /* This matching of child components of v1 and child
                         * components of v2 is invalid, because there are no
                         * isomorphisms of child component k of v1 with child
                         * component comp_perm[k] of v2 */
                        valid = false;
                        break;
                    }
                    score += tmp_max_score;
                    size += tmp_max_size;
                }
                if (!valid)
                    continue;
                /* Check if this comp_perm is the best so far */
                if (score > max_score || (score == max_score
                            && size > max_size)) {
                    max_score = score;
                    max_size = size;
                    max_perms = std::vector<const std::vector<int>* >(1,
                            &comp_perm);
                    max_option_inds = std::vector<std::vector<
                        std::vector<int> > >(1, option_inds);
                } else if (score == max_score && size == max_size) {
                    max_perms.push_back(&comp_perm);
                    max_option_inds.push_back(option_inds);
                }
            }
            final_score += max_score;
            final_size += max_size;
            if (keep_options) {
                /* Update the list of isomorphisms that can attain the best
                 * score and size */
                std::vector<DPNodes> new_options;
                for (unsigned j = 0; j < max_perms.size(); ++j) {
                    std::vector<DPNodes> tmp_options = options;
                    for (unsigned k = 0; k < max_perms[j]->size(); ++k) {
                        if ((*max_perms[j])[k] == -1) continue;
                        options_append(tmp_options, te1[k].first,
                                te2[(*max_perms[j])[k]].first,
                                max_option_inds[j][k]);
                    }
                    new_options.insert(new_options.end(), tmp_options.begin(),
                            tmp_options.end());
                }
                options = new_options;
            }
        }
    }
}

/* Determine the best subgraph isomorphism(s) of mol1 and mol2 that map
 * biconnected component root1 of mol1 to biconnected component root2 of mol2.
 *
 * First, the graphs of biconnected components of mol1 and mol2 (tree1 and tree2
 * respectively) are treated as rooted trees with roots root1 and root2. Then,
 * the following dynamic programming algorithm is applied:
 *
 * For each v1, v2, and iso_ind, best_matches[v1][v2][iso_ind] keeps track of
 * the best isomorphisms of a subtree (in tree1) rooted at v1 with a subtree (in
 * tree2) rooted at v2, such that the isomorphism of components v1 and v2 is the
 * iso_ind'th isomorphism in iso_mat[v1][v2]. best_matches is updated from the
 * components furthest from root1 and root2 down to the roots.
 *
 * When this function returns, for each v1 and v2 that are descendants of root1
 * and root2 respectively,
 * (1) best_matches[v1][v2][iso_ind].score is the highest total score of any
 *     isomorphism matching a subtree of v1 to a subtree of v2, such that the
 *     iso_ind'th isomorphism is used to match v1 and v2.
 * (2) best_matches[v1][v2][iso_ind].size is the largest number of matched atoms
 *     among all isomorphisms yielding the maximal score in (1).
 * (3) best_matches[v1][v2][iso_ind].options is a list of all possible
 *     matchings between children of v1 and children of v2 that yield the
 *     maximal score of (1) and the maximal number of matched atoms of (2). This
 *     is updated only if keep_options is true.
 *
 * The maximal score in (1) and the maximal number of matched atoms for that 
 * score in (2) are returned.
 */
std::pair<double, int> helpers::max_match(const GraphRepr& tree1,
        const TreeEdges& tree_edges1,
        const std::vector<GraphRepr>& components1,
        const std::vector<std::vector<int> > components_idx1, int root1,
        const GraphRepr& tree2, const TreeEdges& tree_edges2,
        const std::vector<GraphRepr>& components2,
        const std::vector<std::vector<int> > components_idx2, int root2,
        const IsoMat& iso_mat, DPMat& best_matches, bool keep_options) {
    if (iso_mat[root1][root2].size() == 0)
        return std::make_pair(0,0);
    /* Treat tree1 and tree2 as rooted with roots root1 and root2.
     *
     * We take a shortcut here and ignore nodes of tree1 with index < root1, as
     * if the best matching involves such a node, it would have been considered
     * in a previous call to max_match. */
    std::vector<std::vector<int> > levels1;
    std::vector<int> level_map1;
    std::vector<std::vector<int> > levels2;
    std::vector<int> level_map2;
    get_levels(tree1, root1, levels1, level_map1, true);
    get_levels(tree2, root2, levels2, level_map2, false);
    /* Dynamic programming algorithm */
    for (int j = levels1.size() - 1; j >= 0; --j) {
        if (j >= int(levels2.size())) continue;
        /* Update best_matches[v1][v2] for all v1 and v2 at distance j from
         * root1 and root2 */
        BOOST_FOREACH(int v1, levels1[j]) {
            BOOST_FOREACH(int v2, levels2[j]) {
                if (iso_mat[v1][v2].size() == 0) continue;
                dp_update(v1, v2, best_matches, tree_edges1[v1], level_map1,
                        tree_edges2[v2], level_map2, iso_mat, keep_options);
            }
        }
    }
    /* Return best score and size over all isomorphisms of root1 and root2 */
    double max_score = -std::numeric_limits<double>::max();
    int max_size = 0;
    BOOST_FOREACH(const DPBest& best, best_matches[root1][root2]) {
        if (best.score > max_score ||
                (best.score == max_score && best.size > max_size)) {
            max_score = best.score;
            max_size = best.size;
        }
    }
    return std::make_pair(max_score, max_size);
}

/* Enumerate all atom maps that attain the best score and most matched atoms
 * for that score, given a best_matches array where best_matches[v1][v2] has
 * been completed using dynamic programming.
 *
 * This function is currently not used; keep it around in case we switch the
 * method later.
 */
#if 0
void recursive_append_matches(const DPMat& best_matches, int v1, int v2, int iso_ind,
        const std::vector<std::vector<int> >& components_idx1,
        const std::vector<int>& atom_idx1,
        const std::vector<std::vector<int> >& components_idx2,
        const std::vector<int>& atom_idx2, const IsoMat& iso_mat,
        std::vector<MatchList>& matches) {
    BOOST_FOREACH(const DPNodes& nodes, best_matches[v1][v2][iso_ind].options) {
        std::vector<std::vector<MatchList> > submatches(nodes.size());
        for (unsigned i = 0; i < nodes.size(); ++i) {
            const DPNode& node = nodes[i];
            recursive_append_matches(best_matches, node.v1,
                        node.v2, node.iso_ind, components_idx1, atom_idx1,
                        components_idx2, atom_idx2, iso_mat, submatches[i]);
        }
        std::vector<int> inds(nodes.size(), 0);
        while (true) {
            MatchList combined;
            std::vector<bool> done;
            for (unsigned i = 0; i < iso_mat[v1][v2][iso_ind].perm.size(); ++i) {
                int a = atom_idx1[components_idx1[v1][i]];
                if (a >= int(done.size()))
                    done.resize(a + 100, false);
                if (done[a]) continue;
                done[a] = true;
                int b = atom_idx2[components_idx2[v2][
                    iso_mat[v1][v2][iso_ind].perm[i]]];
                combined.push_back(std::make_pair(a,b));
            }
            for (unsigned i = 0; i < nodes.size(); ++i) {
                for (unsigned j = 0; j < submatches[i][inds[i]].size(); ++j) {
                    int a = submatches[i][inds[i]][j].first;
                    if (a >= int(done.size()))
                        done.resize(a + 100, false);
                    if (done[a]) continue;
                    done[a] = true;
                    combined.push_back(submatches[i][inds[i]][j]);
                }
            }
            matches.push_back(combined);
            int i = 0;
            while (i < int(submatches.size()) && inds[i] == int(submatches[i].size() - 1)) {
                inds[i] = 0;
                ++i;
            }
            if (i == int(submatches.size())) break;
            ++inds[i];
        }
    }
}
#endif

/* Compute the squared distance between 3-vectors (p1 * rot + trans) and p2 */
double compute_dist(const double* p1, const double* p2, double* rot,
        double* trans) {
    double dx = p2[0]-(rot[0]*p1[0]+rot[3]*p1[1]+rot[6]*p1[2]+trans[0]);
    double dy = p2[1]-(rot[1]*p1[0]+rot[4]*p1[1]+rot[7]*p1[2]+trans[1]);
    double dz = p2[2]-(rot[2]*p1[0]+rot[5]*p1[1]+rot[8]*p1[2]+trans[2]);
    return (dx*dx+dy*dy+dz*dz);
}

/* Compute the best RMSD alignment of a list of matched atoms of mol1 and mol2.
 *
 * If pos1 and pos2 are the Nx3 matrices of positions of atoms of mol1 and mol2,
 * then rot and trans are the 3x3 and 3x1 rotation and translation matrices
 * such that pos2 ~ pos1 * rot + trans. The RMSD value is returned.
 *
 * Note: We roll our own implementation here using the C++ Eigen library, as the
 * SVD implementation in periodicfix does not seem to correctly handle singular
 * matrices.
 */
double compute_align(SystemPtr mol1, SystemPtr mol2,
        const MatchList& match, double* rot, double* trans) {
    /* Get centroids */
    Eigen::Vector3d xm = Eigen::Vector3d::Zero();
    Eigen::Vector3d ym = Eigen::Vector3d::Zero();
    for (unsigned i = 0; i < match.size(); ++i) {
        xm(0) += mol1->atom(match[i].first).x;
        xm(1) += mol1->atom(match[i].first).y;
        xm(2) += mol1->atom(match[i].first).z;
        ym(0) += mol2->atom(match[i].second).x;
        ym(1) += mol2->atom(match[i].second).y;
        ym(2) += mol2->atom(match[i].second).z;
    }
    xm /= match.size();
    ym /= match.size();
    /* Get centered position matrices */
    Eigen::Matrix<double, Eigen::Dynamic, 3> Xo(match.size(), 3);
    Eigen::Matrix<double, Eigen::Dynamic, 3> Yo(match.size(), 3);
    for (unsigned i = 0; i < match.size(); ++i) {
        Xo(i,0) = mol1->atom(match[i].first).x - xm(0);
        Xo(i,1) = mol1->atom(match[i].first).y - xm(1);
        Xo(i,2) = mol1->atom(match[i].first).z - xm(2);
        Yo(i,0) = mol2->atom(match[i].second).x - ym(0);
        Yo(i,1) = mol2->atom(match[i].second).y - ym(1);
        Yo(i,2) = mol2->atom(match[i].second).z - ym(2);
    }
    /* Rotation matrix A maximizes tr(Yo^T*Xo*A) */
    Eigen::Matrix3d prod = Yo.transpose() * Xo;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(prod,
            Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d Ap;
    Ap << 1,1, (svd.matrixU() * svd.matrixV()).determinant();
    Eigen::Matrix3d A = svd.matrixV() * Ap.asDiagonal()
        * svd.matrixU().transpose();
    /* Translation vector b is given by ym - xm*A */
    Eigen::Vector3d b = ym - A.transpose() * xm;
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            rot[i*3+j] = A(i,j);
    for (unsigned i = 0; i < 3; ++i)
        trans[i] = b(i);
    /* RMSD squared is given by (||Xo||+||Yo||-2*tr(Yo^T*Xo*A))/n */
    double xnorm = (Xo.transpose() * Xo).trace();
    double ynorm = (Yo.transpose() * Yo).trace();
    double rmsd_sq = (xnorm+ynorm-2*Ap.dot(svd.singularValues()))/match.size();
    return (rmsd_sq > 0 ? std::sqrt(rmsd_sq) : 0);
}
 
/* Heuristic search algorithm to select, among multiple matches, the match that
 * minimizes RMSD alignment.
 *
 * This function assumes best_matches[v1][v2][iso_ind] has been updated by the
 * dynamic programming algorithm, and it searches for the best atom map that
 * maps component v1 to component v2 by the iso_ind'th isomorphism.
 *
 * match is a list of matching pairs of atom IDs from mol1 and mol2. Returns the
 * minimum RMSD value.
 */
double helpers::select_single(const DPMat& best_matches, int v1, int v2,
        int iso_ind, const std::vector<std::vector<int> >& components_idx1,
        const std::vector<int>& atom_idx1, SystemPtr mol1,
        const std::vector<std::vector<int> >& components_idx2,
        const std::vector<int>& atom_idx2, SystemPtr mol2,
        const IsoMat& iso_mat, MatchList& match) {
    match.clear();

    /* First pass: Find atoms of mol1 for which there is a unique choice of
     * atom in mol2 under any "best score and size" mapping. */
    std::set<int> done_components;
    std::set<int> done_atoms;
    std::deque<const DPBest*> q;
    q.push_back(&best_matches[v1][v2][iso_ind]);
    done_components.insert(v1);
    for (unsigned i = 0; i < iso_mat[v1][v2][iso_ind].perm.size(); ++i) {
        int a = atom_idx1[components_idx1[v1][i]];
        int b = atom_idx2[components_idx2[v2][
            iso_mat[v1][v2][iso_ind].perm[i]]];
        if (done_atoms.find(a) == done_atoms.end()) {
            done_atoms.insert(a);
            match.push_back(std::make_pair(a,b));
        }
    }
    while (q.size() > 0) {
        const DPBest* front = q.front();
        q.pop_front();
        /* Count the number of different mappings to mol2 that each child
         * component of front in mol1 can map to, and vice versa */
        std::map<int, std::set<std::pair<int, int> > > options_per_comp1;
        std::map<int, std::set<std::pair<int, int> > > options_per_comp2;
        BOOST_FOREACH(const DPNodes& option, front->options) {
            BOOST_FOREACH(const DPNode& node, option) {
                options_per_comp1[node.v1].insert(
                        std::make_pair(node.v2, node.iso_ind));
                options_per_comp2[node.v2].insert(
                        std::make_pair(node.v1, node.iso_ind));
            }
        }
        for (std::map<int, std::set<std::pair<int, int> > >::iterator
                iter = options_per_comp1.begin();
                iter != options_per_comp1.end(); ++iter) {
            if (iter->second.size() > 1) continue;
            int v2 = iter->second.begin()->first;
            if (options_per_comp2[v2].size() > 1) continue;
            /* If there is only one possibility, add the matched atoms and
             * continue traversing this subtree */
            int new_v1 = iter->first;
            int new_v2 = iter->second.begin()->first;
            int new_iso = iter->second.begin()->second;
            q.push_back(&best_matches[new_v1][new_v2][new_iso]);
            done_components.insert(new_v1);
            for (unsigned i = 0; i < iso_mat[new_v1][new_v2][
                    new_iso].perm.size(); ++i) {
                int a = atom_idx1[components_idx1[new_v1][i]];
                int b = atom_idx2[components_idx2[new_v2][
                    iso_mat[new_v1][new_v2][new_iso].perm[i]]];
                if (done_atoms.find(a) == done_atoms.end()) {
                    done_atoms.insert(a);
                    match.push_back(std::make_pair(a,b));
                }
            }
        }
    }

    /* Compute the RMSD alignment of the atoms matched in the first pass */
    double rot[9];
    double trans[3];
    compute_align(mol1, mol2, match, rot, trans);

    /* Second pass: At each branch for which a component of mol1 can be mapped
     * to multiple components of mol2, select the mapping that gives the lowest
     * RMSD under the current alignment. */
    q.clear();
    q.push_back(&best_matches[v1][v2][iso_ind]);
    while (q.size() > 0) {
        const DPBest* front = q.front();
        q.pop_front();
        double min_dist = std::numeric_limits<double>::max();
        const DPNodes* best_option = NULL;
        BOOST_FOREACH(const DPNodes& option, front->options) {
            double dist = 0;
            BOOST_FOREACH(const DPNode& node, option) {
                const Isomorphism& iso = iso_mat[node.v1][node.v2][
                    node.iso_ind];
                for (unsigned i = 0; i < iso.perm.size(); ++i) {
                    int a = atom_idx1[components_idx1[node.v1][i]];
                    int b = atom_idx2[components_idx2[node.v2][iso.perm[i]]];
                    double tmp_dist = compute_dist(mol1->atom(a).pos(),
                            mol2->atom(b).pos(), rot, trans);
                    dist += tmp_dist;
                }
            }
            if (dist < min_dist) {
                min_dist = dist;
                best_option = &option;
            }
        }
        if (best_option == NULL)
            MSYS_FAIL("Error -- uninitialized best option");
        BOOST_FOREACH(const DPNode& node, *best_option) {
            q.push_back(&best_matches[node.v1][node.v2][node.iso_ind]);
            if (done_components.find(node.v1) != done_components.end())
                continue;
            done_components.insert(node.v1);
            const Isomorphism& iso = iso_mat[node.v1][node.v2][node.iso_ind];
            for (unsigned i = 0; i < iso.perm.size(); ++i) {
                int a = atom_idx1[components_idx1[node.v1][i]];
                int b = atom_idx2[components_idx2[node.v2][iso.perm[i]]];
                if (done_atoms.find(a) == done_atoms.end()) {
                    done_atoms.insert(a);
                    match.push_back(std::make_pair(a,b));
                }
            }
        }
        /* Update the alignment to include the new mapped atoms */
        compute_align(mol1, mol2, match, rot, trans);
    }
    return compute_align(mol1, mol2, match, rot, trans);
}

void atommatch::AtomMatch(SystemPtr mol1,
        SystemPtr mol2, ScoreFctPtr rep, MatchList& match) {
    /* Initialize graphs from molecules */
    std::vector<int> atom_idx1;
    std::vector<int> atom_idx2;
    GraphRepr g1;
    GraphRepr g2;
    init_graph(mol1, atom_idx1, g1);
    init_graph(mol2, atom_idx2, g2);

    /* Get trees of biconnected components */
    std::vector<GraphRepr> components1;
    std::vector<std::vector<int> > components_idx1;
    GraphRepr tree1;
    TreeEdges tree_edges1;
    std::vector<GraphRepr> components2;
    std::vector<std::vector<int> > components_idx2;
    GraphRepr tree2;
    TreeEdges tree_edges2;
    partition_graph_biconnected(g1, components1, components_idx1, tree1,
            tree_edges1);
    partition_graph_biconnected(g2, components2, components_idx2, tree2,
            tree_edges2);

    /* Get matrix of isomorphisms and scores between all components of g1 and
     * g2 */
    IsoMat iso_mat(components1.size(),
            std::vector<Isomorphisms>(components2.size()));
    for (unsigned i = 0; i < components1.size(); ++i)
        for (unsigned j = 0; j < components2.size(); ++j)
            isomorphisms(components1[i], components_idx1[i], atom_idx1,
                    mol1, components2[j], components_idx2[j], atom_idx2, mol2,
                    rep, iso_mat[i][j]);

    /* Compute the maximal score and maximal size for that score, over all
     * subgraph isomorphisms */
    double max_score = 0;
    int max_size = 0;
    std::vector<std::pair<int, int> > max_roots;
    for (unsigned i = 0; i < components1.size(); ++i) {
        for (unsigned j = 0; j < components2.size(); ++j) {
            DPMat best_matches(tree1.v_to_e.size(),
                    std::vector<std::vector<DPBest> >(tree2.v_to_e.size()));
            /* First time around, do not keep track of isomorphisms */
            std::pair<double, int> output = max_match(tree1, tree_edges1,
                    components1, components_idx1, i, tree2, tree_edges2,
                    components2, components_idx2, j, iso_mat, best_matches,
                    false);
            if (output.first > max_score ||
                    (output.first == max_score && output.second > max_size)) {
                max_score = output.first;
                max_size = output.second;
                max_roots.clear();
                max_roots.push_back(std::make_pair(i,j));
            } else if (output.first == max_score && output.second == max_size)
                max_roots.push_back(std::make_pair(i,j));
        }
    }

    match.clear();
    if (max_score == 0 && max_size == 0)
        return;

    /* Perform a heuristic search for the match with the minimal RMSD among
     * isomorphisms that attain the maximal score and size for that score */
    double min_rmsd = std::numeric_limits<double>::max();
    for (unsigned i = 0; i < max_roots.size(); ++i) {
        DPMat best_matches(tree1.v_to_e.size(),
                std::vector<std::vector<DPBest> >(tree2.v_to_e.size()));
        /* Second time around, keep track of isomorphisms */
        max_match(tree1, tree_edges1, components1, components_idx1,
                max_roots[i].first, tree2, tree_edges2, components2,
                components_idx2, max_roots[i].second, iso_mat, best_matches,
                true);
        for (unsigned j = 0; j < best_matches[max_roots[i].first][
                max_roots[i].second].size(); ++j) {
            if (best_matches[max_roots[i].first][
                    max_roots[i].second][j].score == max_score
                    && best_matches[max_roots[i].first][
                    max_roots[i].second][j].size == max_size) {
                MatchList tmp_match;
                double rmsd = select_single(best_matches, max_roots[i].first,
                        max_roots[i].second, j, components_idx1, atom_idx1,
                        mol1, components_idx2, atom_idx2, mol2, iso_mat,
                        tmp_match);
                if (rmsd < min_rmsd) {
                    min_rmsd = rmsd;
                    match = tmp_match;
                }
            }
        }
    }
}
