#ifndef desres_pfx_graph_hxx
#define desres_pfx_graph_hxx

#include <vector>
#include <algorithm>

namespace desres { namespace msys { namespace pfx {

    // Graph is a helper class for assembling bonds into a graph
    // structure.
    class Graph {
        typedef std::vector<unsigned> IdList;
        std::vector<IdList> g;
        unsigned edgecount;

    public:
        typedef unsigned Id;

        explicit Graph(unsigned nverts)
        : g(nverts), edgecount(0) {}

        // number of vertices
        unsigned nverts() const { return g.size(); }
        
        // number of edges
        unsigned nedges() const { return edgecount; }

        IdList const& atoms(Id i) const { return g.at(i); }

        // Add an edge to the graph between specified vertices.  i==j or
        // i>=size are ignored.  Return true if added or already present;
        // false if unable to add.
        bool add_edge(Id vi, Id vj) {
            if (vi==vj || vi>=g.size() || vj>=g.size()) return false;
            IdList& bi = g[vi];
            IdList& bj = g[vj];
            if (std::find(bj.begin(), bj.end(), vi)==bj.end()) {
                ++edgecount;
                bj.push_back(vi);
                bi.push_back(vj);
            }
            return true;
        }

        // copy_bonds assembles a packed list of bonds in topologically
        // sorted order.  dst should be a back_inserter, or a simple
        // iterator if space is pre-allocated with 2*nedges elements.
        template <typename T>
        void copy_bonds(T dst) const {
            IdList stack;
            std::vector<bool> flags(g.size());
            for (Id i=0, n=g.size(); i<n; i++) {
                stack.push_back(i);
                while (!stack.empty()) {
                    Id ai = stack.back();
                    stack.pop_back();
                    if (flags[ai]) continue;
                    IdList const& b = g[ai];
                    stack.insert(stack.end(), b.begin(), b.end());
                    for (unsigned j=0, m=b.size(); j<m; j++) {
                        Id aj=b[j];
                        if (!flags[aj]) {
                            *dst++ = ai;
                            *dst++ = aj;
                        }
                    }
                    flags[ai]=true;
                }
            }
        }

        // copy_components assembles packed lists of vertices grouped by 
        // connected component, and the sizes of the corresponding components.
        template <typename T>
        void copy_components(T verts, T sizes) const {
            IdList stack;
            std::vector<bool> flags(g.size());
            for (unsigned i=0, n=g.size(); i<n; i++) {
                stack.push_back(i);
                Id size=0;
                while (!stack.empty()) {
                    Id ai = stack.back();
                    stack.pop_back();
                    if (flags[ai]) continue;
                    IdList const& b = g[ai];
                    stack.insert(stack.end(), b.begin(), b.end());
                    *verts++ = ai;
                    ++size;
                    flags[ai]=true;
                }
                if (size!=0) *sizes++ = size;
            }
        }
    };

}}}

#endif
