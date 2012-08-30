#ifndef desres_msys_analyze_aromatic_hxx
#define desres_msys_analyze_aromatic_hxx

#include "../system.hxx"

namespace desres { namespace msys {
    struct AromaticAtomClassification {
        enum {
            X_TYPE=0,
            Y_TYPE,
            YEXT_TYPE,
            Z_TYPE,
            INVALID
        };
    };
    struct AromaticRingClassification {
        enum {
            ANTIAROMATIC=0,
            NONAROMATIC,
            AROMATIC
        };
    };

    /* Balaban, Chem Rev, 2004, 104, 2777-2812 
       Implementation of *simple* aromaticity detection
       Given ring atom types {X[2e], Y[1e], Z[0e]}, a ring of size M should satisfy:
       X+Y+Z=M
       X+0.5Y=2n+1  -> aromatic
       X+0.5Y=2n    -> antiaromatic
       The following additional restrictions are *NOT* considered here:
       1) sum of aromaticity constants should be -200< k[A] <200
       2) adjacent atoms of same type (for X and Z) destabilize the ring 
    */
    unsigned classify_ring_atom(unsigned nb, unsigned a0, unsigned b0, unsigned b1, unsigned be);
    unsigned classify_ring_aromaticity(unsigned nX, unsigned nY, unsigned nYe, unsigned nZ);
    unsigned classify_ring_aromaticity(SystemPtr mol, IdList const& atoms);

    double ring_planarity_descriptor(SystemPtr mol, IdList const& aids);

}}

#endif

