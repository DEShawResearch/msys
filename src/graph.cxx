#include "graph.hxx"
#include <stack>
#include <queue>
#include <sstream> // for ostringstream
#include <iostream> // for debugging

using namespace desres::msys;

std::string Graph::hash(SystemPtr sys, IdList const& atoms){
    /* Compute formula of a set of atoms. */
    static const int max_atomic_number=128;
    std::vector<Id> countmap(max_atomic_number,0);
    int biggest=0;
    Id bcount=0;
    for (Id id : atoms){
        int anum=sys->atom(id).atomic_number;
        if (anum < 1) continue;
        ++countmap.at(anum);
        const IdList& bonded = sys->bondedAtoms(id);
        for (unsigned i = 0; i < bonded.size(); ++i) {
            if (sys->atom(bonded[i]).atomic_number != 0)
                ++bcount;
        }
        if (anum>biggest) biggest=anum;
    }
    std::ostringstream ss("");
    // bcount = 2 * (# internal bonds) + (# external bonds)
    ss << bcount;
    
    for (int i=1; i<=biggest; ++i){
        if(countmap[i]){
            ss << " " << i << " " << countmap[i];
        }
    }
    return ss.str();
}

GraphPtr Graph::create(SystemPtr sys, const IdList& atoms) {
    return GraphPtr(new Graph(sys, atoms));
}

Graph::Graph(SystemPtr sys, const IdList& atoms) {
    typedef std::map<std::pair<int,Id>, Id> attrHash;
    typedef std::multimap<Id,Id> countMap;

    _sys = sys;
    attrHash hash_to_idx;
    countMap count_to_idx;
    msys::MultiIdList freq_partition;
    /* Initialize frequency of atom occurance */
    for (Id id : atoms){
        atom_t const& atm = sys->atom(id);
        if (atm.atomic_number < 1) continue;    
        std::pair<int,Id> key(atm.atomic_number, sys->bondCountForAtom(id));
        attrHash::iterator ihash=hash_to_idx.lower_bound(key);
        if (ihash==hash_to_idx.end() || hash_to_idx.key_comp()(key, ihash->first) ){
           ihash=hash_to_idx.insert(ihash, attrHash::value_type(key, freq_partition.size()));
           freq_partition.push_back(msys::IdList());
        }
        freq_partition[ihash->second].push_back(id);
    }
    /* count em up */
    for(Id idx=0; idx<freq_partition.size();++idx){
        count_to_idx.insert(countMap::value_type(freq_partition[idx].size(),idx));
    }

    int prefix_deg = 0;
    /* mapping from id back to index within the node list */
    std::map<Id, int> id_map;

    /* first pass: construct the offsets into the nbr table */
    for (auto const& entry : count_to_idx){
        for (Id id : freq_partition[entry.second]){
            /* Already checked above 
             * atom_t const& atm = sys->atom(id);
             * if (atm.atomic_number < 1) continue; 
             */
            id_map[id] = _nodes.size();
            _ids.push_back(id);
            _nodes.push_back(Node());
            Node& node = _nodes.back();
            node.nnbr = prefix_deg;
            prefix_deg += sys->bondCountForAtom(id);
        }
    }
    _nbrs.resize(prefix_deg);
    /* second pass: update the nbr pointer in each node to point into
     * the nbr table, find the neighbors of each node that are within
     * the atom set, and compute the node attribute as the atomic number
     * combined with the number of non-pseudo bonds */
    int nidx=0;
    for (auto const& entry : count_to_idx){
        for (Id id : freq_partition[entry.second]){
            atom_t const& atm = sys->atom(id);
            /* Already checked above
             * if (atm.atomic_number < 1) continue;
             */
            Node& node = _nodes[nidx++];
            node.nbr = &_nbrs[node.nnbr];
            node.nnbr = 0;
            int degree = 0;
            IdList const& bonded = sys->bondedAtoms(id);
            for (Id other : bonded) {
                if (sys->atom(other).atomic_number == 0) continue; // Pseudo atom
                ++degree;
                if (sys->atom(other).atomic_number == -1) continue; // External atom
                std::map<Id,int>::const_iterator it = id_map.find(other);
                if (it != id_map.end())
                    node.nbr[node.nnbr++] = it->second;
            }
            node.attr = (atm.atomic_number << 16) | degree;
        }
    }
}

void Graph::setNodeAttributes(std::vector<int> const& newattr){
    assert(newattr.size()==_nodes.size());
    for(size_t idx=0; idx<newattr.size();++idx){
        _nodes[idx].attr=newattr[idx];
    }
}


bool Graph::match(const GraphPtr other, std::vector<IdPair>& perm) const {
    perm.clear();
    if (size() != other->size())
        return false;
    if (size() == 0) {
        perm.clear();
        return true;
    }
    for (unsigned i = 0; i < size(); ++i) {
        if (match_common(other, 0, i, perm))
            return true;
    }
    return false;
}

bool Graph::match(const GraphPtr other, msys::Id this_root, msys::Id other_root,
        std::vector<IdPair>& perm) const {
    perm.clear();
    if (size() != other->size())
        return false;
    int this_node = -1;
    int other_node = -1;
    for (unsigned i = 0; i < size(); ++i) {
        if (_ids[i] == this_root)
            this_node = i;
        if (other->_ids[i] == other_root)
            other_node = i;
    }
    if (this_node == -1 || other_node == -1)
        MSYS_FAIL("Invalid root node specifications for graph isomorphism");
    return match_common(other, this_node, other_node, perm);
}


/* Match whether (1) node g of Graph G matches node h of Graph h in attribute,
 * (2) already matched neighbors of g correspond with neighbors of h, and
 * (3) already matched neighbors of h correspond with neighbors of g */
bool Graph::match_node(const Node& g, const Node& h, const std::vector<int>&
       GtoH, const std::vector<int>& HtoG, bool include_nnbr) {
   if (g.attr != h.attr)
       return false;
   if (include_nnbr && g.nnbr != h.nnbr)
       return false;
   for (int i = 0; i < g.nnbr; ++i) {
       int hnbr = GtoH[g.nbr[i]];
       if (hnbr == -1)
           continue;
       /* Check that matched neighbors of g map to neighbors of h */
       bool is_nbr = false;
       for (int j = 0; j < h.nnbr; ++j) {
           if (h.nbr[j] == hnbr)
               is_nbr = true;
       }
       if (!is_nbr)
           return false;
   }
   for (int j = 0; j < h.nnbr; ++j) {
       int gnbr = HtoG[h.nbr[j]];
       if (gnbr == -1)
           continue;
       /* Check that matched neighbors of h map to neighbors of g */
       bool is_nbr = false;
       for (int i = 0; i < g.nnbr; ++i) {
           if (g.nbr[i] == gnbr)
               is_nbr = true;
       }
       if (!is_nbr)
           return false;
   }
   return true;
}

struct Quadruple {
    Quadruple(int _g, int _g_nbr_ind, int _h, int _h_nbr_ind) :
        g(_g), g_nbr_ind(_g_nbr_ind), h(_h), h_nbr_ind(_h_nbr_ind) { }
    int g;
    int g_nbr_ind;
    int h;
    int h_nbr_ind;

    /* For debugging */
    friend std::ostream& operator<<(std::ostream& stream, const Quadruple& q) {
        stream << "(" << q.g << "," << q.g_nbr_ind << "," << q.h
            << "," << q.h_nbr_ind << ")";
        return stream;
    }
};

/* Match this graph to other graph, requiring that this_root matches
 * other_root */
bool Graph::match_common(const GraphPtr other, int this_root, int other_root,
        std::vector<IdPair>& perm) const {
    if (size() != other->size())
        return false;
    if (_nodes[this_root].nnbr != other->_nodes[other_root].nnbr)
        return false;
    if (_nodes[this_root].attr != other->_nodes[other_root].attr)
        return false;
    /* Root nodes match */

    /* Mappings between node in G (this graph) and node in H (other graph),
     * with default (unmatched) set to -1 */
    unsigned n = size();
    std::vector<int> GtoH(n,-1);
    std::vector<int> HtoG(n,-1);
    GtoH[this_root] = other_root;
    HtoG[other_root] = this_root;
    if (n == 1) {
        assert(this_root == 0 && other_root == 0);
        perm.push_back(std::make_pair(_ids[0], other->_ids[0]));
        return true;
    }

    /* Stacks of already matched nodes in G and H */
    std::stack<int> matched_G;
    std::stack<int> matched_H;
    matched_G.push(this_root);
    matched_H.push(other_root);

    /* Stack of (G_node, G_nbr_index, H_node, H_nbr_index) where G_node has
     * been matched with H_node and the G_nbr_index'th neighbor of G_node has
     * been matched with the H_nbr_index'th neighbor of H_node */
    std::stack<Quadruple> progress; 

    Quadruple q(this_root, 0, other_root, -1);
    while (true) {
        /* Find match for neighbor q.g_nbr_ind of q.g */
        int gnbr = _nodes[q.g].nbr[q.g_nbr_ind];
        int hnbr = -1;
        for (int j = q.h_nbr_ind + 1; j < other->_nodes[q.h].nnbr; ++j) {
            /* Does neighbor j of q.h match? */
            if (HtoG[other->_nodes[q.h].nbr[j]] != -1)
                continue; /* No, neighbor j of q.h was already matched */
            if (!match_node(_nodes[gnbr],
                        other->_nodes[other->_nodes[q.h].nbr[j]], GtoH, HtoG))
                continue; /* No, neighbor j of q.h does not match */
            /* Yes, neighbor j of q.h matches */
            q.h_nbr_ind = j;
            hnbr = other->_nodes[q.h].nbr[j];
            break;
        }
        if (hnbr != -1) {
            /* Match was found; node gnbr matches node hnbr */
            progress.push(q);
            GtoH[gnbr] = hnbr;
            HtoG[hnbr] = gnbr;
            matched_G.push(gnbr);
            matched_H.push(hnbr);
            if (matched_G.size() == n) {
                /* All nodes matched */
                for (unsigned i = 0; i < n; ++i)
                    perm.push_back(std::make_pair(_ids[i],
                                other->_ids[GtoH[i]]));
                return true;
            }

            q = Quadruple(gnbr, -1, hnbr, -1);
            bool found_unmatched;
            do {
                /* Find next unmatched neighbor of q.g */
                found_unmatched = false;
                for (int i = q.g_nbr_ind + 1; i < _nodes[q.g].nnbr; ++i) {
                    if (GtoH[_nodes[q.g].nbr[i]] == -1) {
                        q.g_nbr_ind = i;
                        found_unmatched = true;
                        break;
                    }
                }
                if (!found_unmatched) {
                    /* All neighbors of q.g have been matched. We must back up
                     * the progress stack to find the next unmatched neighbor */
                    if (progress.size() == 0) {
                        /* Somehow there are no unmatched neighbors remaining in
                         * G, but not all the nodes of G have been matched... */
                        MSYS_FAIL("Graph isomorphism error---check that graph"
                                "is connected");
                    }
                    q = progress.top();
                    progress.pop();
                }
            } while (!found_unmatched);
            /* In the next iteration of main loop, we search for a match for
             * q.g's q.g_nbr_ind'th neighbor among all neighbors of q.h */
            q.h_nbr_ind = -1;
        } else {
            /* There is no match for gnbr; we must undo the most recent
             * (g, g_nbr_ind, h, h_nbr_ind) in the progress stack and find the
             * next match for g_nbr_ind */ 
            if (progress.size() == 0) {
                /* No more undos possible---no matches found */
                return false;
            }

            Quadruple last = progress.top();
            /* Undo all matches up to and including the match of
             * last.g_nbr_ind with last.h_nbr_ind */
            int gnbr = _nodes[last.g].nbr[last.g_nbr_ind];
            int hnbr = other->_nodes[last.h].nbr[last.h_nbr_ind];
            while (matched_G.top() != gnbr) {
                GtoH[matched_G.top()] = -1;
                HtoG[matched_H.top()] = -1;
                matched_G.pop();
                matched_H.pop();
            }
            assert(matched_G.top() == gnbr && matched_H.top() == hnbr);
            GtoH[gnbr] = -1;
            HtoG[hnbr] = -1;
            matched_G.pop();
            matched_H.pop();

            /* In the next iteration of the loop, we search for a new match for
             * q.g's q.g_nbr_ind'th neighbor among neighbors of q.h starting
             * with q.h_nbr_ind + 1 */
            q = progress.top();
            progress.pop();
        }
    }
}

/* In principle a small modification of matchAll can implement match and
 * matchCommon above, but this matchAll algorithm is slower than the algorithm
 * of matchCommon */
unsigned Graph::matchAll(const GraphPtr other, std::vector<MatchList>&
        perms, bool substructure) const {
    perms.clear();
    if (size() != other->size() && !substructure)
        return 0;
    if (size() == 0) {
        perms.push_back(MatchList());
        return 1;
    }
    
    /* Nodes of this graph G in BFS order from node 0 */
    std::vector<bool> visited(size(), false);
    std::queue<int> Q;
    std::vector<int> nodes_bfs;
    nodes_bfs.reserve(size());
    Q.push(0);
    nodes_bfs.push_back(0);
    visited[0] = true;
    while (Q.size() > 0) {
        int front = Q.front();
        Q.pop();
        for (int i = 0; i < _nodes[front].nnbr; ++i) {
            if (!visited[_nodes[front].nbr[i]]) {
                Q.push(_nodes[front].nbr[i]);
                nodes_bfs.push_back(_nodes[front].nbr[i]);
                visited[_nodes[front].nbr[i]] = true;
            }
        }
    }
    if (nodes_bfs.size() != size())
        MSYS_FAIL("Graph isomorphism error---check that graph is connected");

    /* Maintain the following: For all but top element of matches, matches[i] of
     * H matches nodes_bfs[i] of G. Top element of matches indicates the node in
     * H we are currently trying to match to nodes_bfs[i] of G. */
    std::stack<unsigned> matches;
    matches.push(0);
    std::vector<int> GtoH(size(), -1);
    std::vector<int> HtoG(other->size(), -1);
    while (matches.size() > 0) {
        int g = nodes_bfs[matches.size() - 1];
        int h = matches.top();
        /* If we are matching a substructure, we do not require the number of
         * internal bonds of g and h (nnbr) to match. Otherwise we do. */
        if (HtoG[h] == -1 && match_node(_nodes[g],
                    other->_nodes[h], GtoH, HtoG, !substructure)) {
            /* Matches; try to match next node of G in nodes_bfs */
            GtoH[g] = h;
            HtoG[h] = g;
            matches.push(0);
        } else {
            /* Does not match; increment top element of matches */
            matches.pop();
            matches.push(h+1);
        }
        if (matches.size() == size() + 1) {
            /* Have matched all nodes in G */
            MatchList perm;
            for (unsigned i = 0; i < size(); ++i)
                perm.push_back(std::make_pair(_ids[i],
                            other->_ids[GtoH[i]]));
            perms.push_back(perm);

            /* Top element is the placeholder 0; remove it */
            matches.pop();

            /* Undo top match and increment */
            h = matches.top();
            GtoH[HtoG[h]] = -1;
            HtoG[h] = -1;
            matches.pop();
            matches.push(h+1);
        }
        while (matches.top() == other->size()) {
            /* Top element is size of H; remove it to return to previous node of
             * G in nodes_bfs */
            matches.pop();
            if (matches.size() == 0)
                break;
            /* Undo top match and increment */
            h = matches.top();
            GtoH[HtoG[h]] = -1;
            HtoG[h] = -1;
            matches.pop();
            matches.push(h+1);
        }
    }
    return perms.size();
}
