#ifndef desres_msys_graph_hxx
#define desres_msys_graph_hxx

#include "system.hxx"
#include <vector>
#include <map>
#include <boost/noncopyable.hpp>

namespace desres { namespace msys {

    /* Checks graph isomorphism for two sets of atoms */ 
    class Graph : boost::noncopyable {

    private:
        struct Node {
            int attr;   /* two bytes for atomic number, two for degree */
            int nnbr;   /* bonds within the template */
            int* nbr;   /* 0-based index of neighbors */
        };

        std::vector<Node> _nodes; /* nodes */
        std::vector<int> _nbrs;  /* storage for neighbors */
        IdList _ids; /* keep track of original IDs */

        /* Helper functions for isomorphism match */
        static bool match_node(const Node& g, const Node& h,
                const std::vector<int>& GtoH, const std::vector<int>& HtoG);
        bool match_common(const boost::shared_ptr<Graph> other, int this_root,
                int other_root, std::vector<IdPair>& matches) const;

        /* Constructor is private; must use create() function */
        Graph(SystemPtr sys, const IdList& atoms);

    public:

        typedef std::vector<IdPair> MatchList;

        /* string hash of attributes in the nodes of the graph */
        static std::string hash(SystemPtr sys, IdList const& atoms);

        /* construct an isomorphism topology using the given atoms.
         * Atoms outside the atom set count towards the degree of
         * the internal atoms but are also not part of the graph.
         * Atoms with atomic_number 0 are ignored. */
        static boost::shared_ptr<Graph> create(SystemPtr sys,
                const IdList& atoms);

        /* number of particles used in the isomorphism check */
        unsigned size() const { return _nodes.size(); }

        /* On success, return true and store a list of (this_id, other_id)
         * matched pairs in matches. Stored IDs are the original atom IDs
         * in the two systems, and external atoms (with atomic number 0) are not
         * included. */
        bool match(const boost::shared_ptr<Graph> other,
                MatchList& matches) const;

        /* Match under the constraint that atom this_root of this graph must
         * match atom other_root of other graph */
        bool match(const boost::shared_ptr<Graph> other, Id this_root,
                Id other_root, MatchList& matches) const;

        /* Find all permutations that match, returning the number of matches
         * and storing list of all matched permutations. */
        unsigned matchAll(const boost::shared_ptr<Graph> other,
                std::vector<MatchList>& matches) const;
        
        /* Ordering of atoms in the graph structure... Used in conjuntion with 
         * setNodeAttributes */
        IdList const& atoms() const {return _ids;}

        /* Overides the default node attribute (int) with user supplied value.
         * Must specify attributes for all atoms in the same order as returned 
         * from atoms() */      
        void setNodeAttributes(std::vector<int> const& newattr);
        
    };
    typedef boost::shared_ptr<Graph> GraphPtr;

}}

#endif
