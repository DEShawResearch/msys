#ifndef desres_msys_aromatic_hxx
#define desres_msys_aromatic_hxx

#include "system.hxx"

namespace desres { namespace msys {

    /* Balaban, Chem Rev, 2004, 104, 2777-2812 
       Implementation of *simple* aromaticity detection
       Given ring atom types {X[2e], Y[1e], Z[0e]}, a ring of size M 
       should satisfy:

       X+Y+Z=M
       X+0.5Y=2n+1  -> aromatic
       X+0.5Y=2n    -> antiaromatic

       The following additional restrictions are *NOT* considered here:
       1) sum of aromaticity constants should be -200< k[A] <200
       2) adjacent atoms of same type (for X and Z) destabilize the ring 
    */

    struct AromaticAtom {
        enum Type {
            X_TYPE=0,
            Y_TYPE,
            YEXT_TYPE,
            Z_TYPE,
            INVALID
        };

        /* FIXME: I have no idea what these inputs mean */
        static Type Classify(int nb, int a0, int b0, int b1, int be);
    };

    struct AromaticRing {
        enum Type {
            ANTIAROMATIC=0,
            NONAROMATIC,
            AROMATIC
        };

        /* Classify according to aromatic atom types */
        static Type Classify(int nX, int nY, int nYe, int nZ);
    };



    /* Find the number of atoms posessing the given AromaticAtomClassification
     * in the given ring.  Return true if all atoms were able to be classified,
     * false if not. */
    bool ClassifyAromaticAtoms(SystemPtr mol, IdList const& atoms,
                               int& nx, int& ny, int& nyext, int& nz);

    /* Get the aromatic ring type of the given atoms. */
    AromaticRing::Type ClassifyAromaticRing(SystemPtr mol, IdList const& atms);

    /* Compute ring planarity. */
    double ComputeRingPlanarity(SystemPtr mol, IdList const& atoms);

    /* A an object which caches structure perception data.
     * For developer use only: Unstable and subject to change without notice.
     */
    class annotation_t : boost::noncopyable {

        SystemPtr _sys;

    public:
        explicit annotation_t(SystemPtr sys);

        /* automatic conversion */
        operator SystemPtr() const { return _sys; }
        System* operator->() const { return _sys.get(); }

        typedef boost::shared_ptr<std::vector<int> > IntListPtr;

        struct atom_data_t {
            unsigned char   aromatic;
            unsigned char   hcount;
            unsigned char   valence;
            unsigned char   degree;
            IntListPtr      ring_sizes;
            unsigned char   ring_bonds;
        
            atom_data_t() 
            : aromatic(), hcount(), valence(), degree(), ring_bonds() {}
        
            int ring_count() const { 
                return ring_sizes ? ring_sizes->size() : 0; 
            }
            bool in_ring(int size = -1) const {
                if (size==-1) return ring_count();
                IntListPtr p = ring_sizes;
                return p && std::find(p->begin(), p->end(), size)!=p->end();
            }
        };
        
        struct bond_data_t {
            bool    aromatic;
            bool    in_ring;
        
            bond_data_t() : aromatic(), in_ring() {}
        };

        std::vector<atom_data_t> atoms;
        std::vector<bond_data_t> bonds;
    
        int hybridization(int ai) const;
    };
}}

#endif

